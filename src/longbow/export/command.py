from cmath import inf
import logging
import math
import sys
import itertools
import re
import time
import os

import click
import click_log
import tqdm

import numpy as np
import pysam
import multiprocessing as mp

import gzip
from construct import *

import longbow.utils.constants
from ..utils import bam_utils
from ..utils.bam_utils import SegmentInfo
from ..utils import model as LongbowModel
from ..utils.model import LibraryModel


logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("export")
click_log.basic_config(logger)


@click.command(name=logger.name)
@click_log.simple_verbosity_option(logger)
@click.option(
    "-p",
    "--pbi",
    required=False,
    type=click.Path(),
    help="BAM .pbi index file",
)
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
    "--output-tag",
    default="-",
    type=click.Path(exists=False),
    help="model name output  [default: stdout]",
)
@click.option(
    "-c",
    "--chunk",
    type=str,
    default="",
    required=False,
    help="Process a single chunk of data (e.g. specify '2/4' to process the second of four equally-sized "
         "chunks across the dataset)"
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    default=False,
    show_default=True,
    help="Force overwrite of the output files if they exist."
)
@click.option(
    "-b",
    "--barcode-tag",
    required=True,
    type=str,
    help="The barcode tag to export"
)
@click.option(
    "-e",
    "--expand-tag",
    type=int,
    default=0,
    help="The barcode tag to export"
)
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(pbi, threads, output_tag, chunk, force, barcode_tag, expand_tag, input_bam):
    """Export a tag from an annotated Longbow BAM for use (e.g. correction) with an external tool."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    # Check to see if the output files exist:
    bam_utils.check_for_preexisting_files(output_tag, exist_ok=force)

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

    pbi = f"{input_bam.name}.pbi" if pbi is None else pbi
    read_count = None
    read_num = 0
    start_offset = 0
    end_offset = math.inf

    if not os.path.exists(pbi) and chunk is not "":
        raise ValueError(f"Chunking specified but pbi file '{pbi}' not found")

    if os.path.exists(pbi):
        if chunk is not "":
            (chunk, num_chunks) = re.split("/", chunk)
            chunk = int(chunk)
            num_chunks = int(num_chunks)

            # Decode PacBio .pbi file and determine the shard offsets.
            offsets, zmw_counts, read_count, read_counts_per_chunk, read_nums = bam_utils.compute_shard_offsets(pbi, num_chunks)

            start_offset = offsets[chunk - 1]
            end_offset = offsets[chunk] if chunk < len(offsets) else offsets[chunk - 1]
            read_count = read_counts_per_chunk[chunk - 1] if chunk < len(offsets) else 0
            read_num = read_nums[chunk - 1] if chunk < len(offsets) else 0

            logger.info("Exporting tags in %d reads from chunk %d/%d (reads %d-%d)", read_count, chunk, num_chunks, read_num, read_num + read_count - 1)
        else:
            read_count = bam_utils.load_read_count(pbi)
            logger.info("Exporting tags in %d reads", read_count)

    # Create queues for data:
    queue_size = threads * 2 if threads < 10 else 20
    manager = mp.Manager()
    input_data_queue = manager.Queue(maxsize=queue_size)
    results = manager.Queue()

    # Start worker sub-processes:
    worker_pool = []

    for i in range(threads):
        p = mp.Process(
            target=_worker_export_barcode_tag_fn, args=(input_data_queue, results, i, barcode_tag, expand_tag)
        )
        p.start()
        worker_pool.append(p)

    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(
        input_bam if start_offset == 0 else input_bam.name, "rb", check_sq=False, require_index=False
    ) as bam_file:
        # If we're chunking, advance to the specified virtual file offset.
        if start_offset > 0:
            bam_file.seek(start_offset)

        # Start output worker:
        res = manager.dict()
        output_worker = mp.Process(
            target=_collect_thread_fn,
            args=(results, output_tag, not sys.stdin.isatty(), res, read_count)
        )
        output_worker.start()

        # Add in a sentinel value at the end of the queue - one for each subprocess - so we guarantee
        # that all subprocesses will exit:
        reads_seen = 0
        iter_data = itertools.chain(bam_file, (None,) * threads)
        for r in iter_data:
            # We have to adjust for our sentinel value if we've got to it:
            if r is not None:
                r = r.to_string()
            input_data_queue.put(r)

            if start_offset > 0:
                if bam_file.tell() >= end_offset or r is None:
                    [input_data_queue.put(None) for _ in range(threads)]
                    break

            # Increment reads seen
            reads_seen += 1

        logger.debug("Finished reading data and sending it to sub-processes.")
        logger.debug("Waiting for sub-processes to finish...")

        # Wait for our input jobs to finish:
        for p in worker_pool:
            p.join()

        logger.debug("All workers stopped.")
        logger.debug("Terminating output process.")

        # Now that our input processes are done, we can add our exit sentinel onto the output queue and
        # wait for that process to end:
        results.put(None)
        output_worker.join()

    # with open(output_model if output_model != "-" else "/dev/stdout", "w") as wm:
    #     wm.write(f'{best_model}\n')

    et = time.time()
    logger.info(f"Done. Elapsed time: {et - t_start:2.2f}s. "
                f"Overall processing rate: {reads_seen/(et - t_start):2.2f} reads/s.")


def _collect_thread_fn(out_queue, out_tag_file_name, disable_pbar, res, read_count):
    """Thread / process fn to write out all our data."""

    with open(out_tag_file_name if out_tag_file_name != "-" else "/dev/stdout", "w") as wm, tqdm.tqdm(
            desc="Progress",
            unit=" read",
            colour="green",
            file=sys.stderr,
            disable=disable_pbar,
            total=read_count
        ) as pbar:

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
            tag = raw_data

            wm.write(f'{tag}\n')

            pbar.update(1)


def _worker_export_barcode_tag_fn(in_queue, out_queue, worker_num, barcode_tag, expand_tag):
    """Function to run in each subthread / subprocess.
    Export selected tag from read."""

    num_reads_processed = 0

    while True:
        # Wait until we get some data.
        # Note: Because we have a sentinel value None inserted at the end of the input data for each
        #       subprocess, we don't have to add a timeout - we're guaranteed each process will always have
        #       at least one element.
        raw_data = in_queue.get()

        # Check for exit sentinel:
        if raw_data is None:
            break
        # Should really never be None, but just in case:
        elif raw_data is None:
            continue

        # Unpack our data here:
        read = raw_data
        read = pysam.AlignedSegment.fromstring(
            read, pysam.AlignmentHeader.from_dict(dict())
        )

        # Process and place our data on the output queue:
        if read.has_tag(barcode_tag):
            tag_value = read.get_tag(barcode_tag)
            out_queue.put(tag_value)

        num_reads_processed += 1

    logger.debug(f"Worker %d: Num reads processed: %d", worker_num, num_reads_processed)
