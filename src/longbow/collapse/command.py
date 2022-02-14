import logging
import math
import sys
import itertools
import re
import time
import os
import io
from collections import defaultdict

import click
import click_log
import tqdm

import ssw

import pysam
import multiprocessing as mp
import subprocess
import tempfile

import gzip
from construct import *

import longbow.utils.constants
from ..utils import bam_utils, barcode_utils
from ..utils.bam_utils import SegmentInfo
from ..utils import model as LongbowModel
from ..utils.model import LibraryModel


logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("collapse")
click_log.basic_config(logger)


@click.command(name=logger.name)
@click_log.simple_verbosity_option(logger)
@click.option(
    "-t",
    "--threads",
    type=int,
    default=mp.cpu_count() - 1,
    show_default=True,
    help="number of threads to use (0 for all)",
)
@click.option(
    "-o",
    "--output-bam",
    default="-",
    type=click.Path(exists=False),
    help="annotated bam output  [default: stdout]",
)
@click.option(
    "-m",
    "--model",
    help="The model to use for annotation.  If not specified, it will be autodetected from "
         "the BAM header.  If the given value is a pre-configured model name, then that "
         "model will be used.  Otherwise, the given value will be treated as a file name "
         "and Longbow will attempt to read in the file and create a LibraryModel from it.  "
         "Longbow will assume the contents are the configuration of a LibraryModel as per "
         "LibraryModel.to_json()."
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    default=False,
    show_default=True,
    help="Force overwrite of the output files if they exist."
)
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(threads, output_bam, model, force, input_bam):
    """Cluster UMIs and sequences."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    # Check to see if the output files exist:
    bam_utils.check_for_preexisting_files(output_bam, exist_ok=force)

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

    # Load barcode allow list
    bc_corrected = _cluster_umis_and_sequences(io.open(input_bam.name, "rb"), threads=threads)

    # Configure process manager:
    # NOTE: We're using processes to overcome the Global Interpreter Lock.
    manager = mp.Manager()
    process_input_data_queue = manager.Queue(threads)
    results = manager.Queue()

    # Start worker sub-processes:
    worker_process_pool = []
    for _ in range(threads):
        p = mp.Process(
            target=_correct_barcode_fn, args=(process_input_data_queue, results, barcode_tag, corrected_tag, bc_corrected)
        )
        p.start()
        worker_process_pool.append(p)

    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(
        input_bam, "rb", check_sq=False, require_index=False
    ) as bam_file, tqdm.tqdm(
        desc="Progress",
        unit=" read",
        colour="green",
        file=sys.stderr,
        leave=False,
        disable=not sys.stdin.isatty(),
    ) as pbar:

        # Get our model:
        if model is None:
            lb_model = LibraryModel.from_json_obj(bam_utils.get_model_from_bam_header(bam_file.header))
        elif model is not None and LibraryModel.has_prebuilt_model(model):
            lb_model = LibraryModel.build_pre_configured_model(model)
        else:
            lb_model = LibraryModel.from_json_file(model)

        logger.info(f"Using %s: %s", lb_model.name, lb_model.description)

        out_header = bam_utils.create_bam_header_with_program_group(logger.name, bam_file.header, models=[lb_model])

        # Start output worker:
        res = manager.dict({"num_reads_corrected": 0, "num_reads": 0})
        output_worker = mp.Process(
            target=_write_thread_fn,
            args=(
                results,
                out_header,
                output_bam,
                pbar,
                res
            ),
        )
        output_worker.start()

        # Add in a sentinel value at the end of the queue - one for each subprocess - so we guarantee
        # that all subprocesses will exit:
        iter_data = itertools.chain(bam_file, (None,) * threads)
        for r in iter_data:
            if r is not None:
                process_input_data_queue.put(r.to_string())
            else:
                process_input_data_queue.put(r)

        # Wait for our input jobs to finish:
        for p in worker_process_pool:
            p.join()

        # Now that our input processes are done, we can add our exit sentinel onto the output queue and
        # wait for that process to end:
        results.put(None)
        output_worker.join()

    logger.info(f"Corrected tags in {res['num_reads_corrected']} reads of {res['num_reads']} total ({100.0*res['num_reads_corrected']/res['num_reads']:.2f}%).")
    
    et = time.time()
    logger.info(f"Done. Elapsed time: {et - t_start:2.2f}s. "
                f"Overall processing rate: {res['num_reads']/(et - t_start):2.2f} reads/s.")


def _write_thread_fn(out_queue, out_bam_header, out_bam_file_name, pbar, res):
    """Thread / process fn to write out all our data."""

    with pysam.AlignmentFile(
        out_bam_file_name, "wb", header=out_bam_header
    ) as out_bam_file:
        while True:
            # Wait for some output data:
            raw_data = out_queue.get()

            # Check for exit sentinel:
            if raw_data is None:
                break
            # Should really never be None, but just in case:
            elif raw_data is None:
                continue

            # Unpack data:
            read, num_segments, num_corrected_segments = raw_data
            read = pysam.AlignedSegment.fromstring(read, out_bam_header)

            # Obligatory log message:
            # logger.debug(
            #     "Path for read %s (%2.2f)%s: %s",
            #     read.query_name,
            #     logp,
            #     " (RC)" if is_rc else "",
            #     segments,
            # )

            # Write our our read:
            out_bam_file.write(read)

            # Increment our counters:
            res["num_reads"] += 1
            if num_corrected_segments > 0:
                res["num_reads_corrected"] += 1

            pbar.update(1)


def _correct_barcode_fn(in_queue, out_queue, barcode_tag, corrected_tag, bc_corrected):
    """Function to run in each subprocess.
    Replace barcode with corrected value."""

    while True:
        # Wait until we get some data.
        # Note: Because we have a sentinel value None inserted at the end of the input data for each
        #       subprocess, we don't have to add a timeout - we're guaranteed each process will always have
        #       at least one element.
        raw_data = in_queue.get()

        # Check for exit sentinel:
        if raw_data is None:
            return

        # Unpack our data here:
        read = pysam.AlignedSegment.fromstring(
            raw_data, pysam.AlignmentHeader.from_dict(dict())
        )

        num_segments = 0
        num_corrected_segments = 0

        if read.has_tag(barcode_tag):
            old_bc = read.get_tag(barcode_tag)
            new_bc = bc_corrected[old_bc] if old_bc in bc_corrected else None

            num_segments += 1
            if new_bc is not None:
                read.set_tag(corrected_tag, new_bc)
                num_corrected_segments += 1

        # Process and place our data on the output queue:
        out_queue.put((read.to_string(), num_segments, num_corrected_segments))


def _cluster_umis_and_sequences(input_bam, pseudocount=1000, threads=1):
    """Extracts UMIs and sequences from input reads and clusters them."""

    # We perform two rounds of clustering on UMI and sequence pairs. This
    # approach follows the algorithm from Guillaume Filion in Starcode.
    # See https://github.com/gui11aume/starcode/blob/master/starcode-umi
    # for furhter details.

    max_len = 1000 # Starcode only clusters sequences up to 1000 bp in length.
    data = [None,]

    uarg = ["starcode", "--seq-id", "-qd2", f"-t{threads}"]
    uproc = subprocess.Popen(uarg, stdout=subprocess.PIPE, stdin=subprocess.PIPE)

    sarg = ["starcode", "--seq-id", "-qd2", f"-t{threads}"]
    sproc = subprocess.Popen(sarg, stdout=subprocess.PIPE, stdin=subprocess.PIPE)

    logger.info("Clustering UMIs...")
    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(
        input_bam, "rb", check_sq=False, require_index=False
    ) as bam_file, tqdm.tqdm(
        desc="Progress",
        unit=" read",
        colour="green",
        file=sys.stderr,
        leave=False,
        disable=not sys.stdin.isatty(),
    ) as pbar:

        for read in bam_file:
            if read.has_tag(longbow.utils.constants.READ_UMI_TAG):
                bc = read.get_tag(longbow.utils.constants.READ_UMI_TAG)
                uproc.stdin.write(bytearray(bc+"\n", 'ascii'))
                sproc.stdin.write(bytearray(read.query_sequence[:max_len]+"\n", 'ascii'))

                data.append((0, read.query_sequence[max_len:]))

            pbar.update(1)

    uo = uproc.communicate()
    so = sproc.communicate()

    # Read UMI clusters
    for line in uo[0].split(b"\n"):
        if line != b'':
            # Parse output.
            centroid, count, ids = line.decode('ascii').rstrip().split("\t")

            # Store UMI cluster ID (same as seq_id of the canonicals).
            for id in ids.split(","):
                data[int(id)] += (centroid,)

    # Read sequence clusters
    cluster_id = 0
    for line in so[0].split(b"\n"):
        if line != b'':
            cluster_id += 1

            # Parse output.
            centroid, count, ids = line.decode('ascii').rstrip().split("\t")

            # Aux vars to recover canonical end.
            canon_end = ""
            seqbuf = {}
            maxcnt = 0

            # Aux vars for final clusters.
            clusters = {}
            for id in ids.split(","):
                # Merge seq id and umi canonical to form final clusters.
                umi = data[int(id)][2]
                if umi not in clusters:
                    clusters[umi] = (id,)
                else:
                    clusters[umi] += (id,)
                    
                seq_end = data[int(id)][1]
                if seq_end not in seqbuf:
                    seqbuf[seq_end] = 1
                else:
                    seqbuf[seq_end] += 1

                if seqbuf[seq_end] > maxcnt:
                    maxcnt = seqbuf[seq_end]
                    canon_end = seq_end

            # Print clusters
            canonical_sequence = centroid + canon_end
            for umi in clusters:
                # Append canonical sequence to UMI cluster.
                clust_out = umi + canonical_sequence

                # Append sequence count.
                clust_out += "\t"+str(len(clusters[umi]))

                # Append sequence ids.
                clust_out += "\t"+",".join(clusters[umi])

                # Write to stdout.
                sys.stdout.write(clust_out.rstrip()+'\n')

    print("")

    return None

