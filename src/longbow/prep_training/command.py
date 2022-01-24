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

import ssw

import pysam
from pysam import FastaFile
import multiprocessing as mp

import gzip
from construct import *

import longbow.utils.constants
from ..utils import bam_utils, barcode_utils
from ..utils.bam_utils import SegmentInfo
from ..utils import model as LongbowModel
from ..utils.model import LibraryModel


logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("prep_training")
click_log.basic_config(logger)


@click.command(name=logger.name)
@click_log.simple_verbosity_option(logger)
@click.option(
    "-s",
    "--sirv-fasta",
    required=True,
    type=click.Path(),
    help="SIRV fasta file",
)
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
# @click.option(
#     "-n",
#     "--training-num-samples",
#     type=int,
#     default=10,
#     show_default=True,
#     help="number of training samples to use",
# )
@click.option(
    "-o",
    "--output-file",
    default="-",
    type=click.Path(exists=False),
    help="training element output file  [default: stdout]",
)
@click.option(
    "-m",
    "--model",
    type=str,
    default=longbow.utils.constants.DEFAULT_MODEL,
    show_default=True,
    help="The model to use for annotation.  If the given value is a pre-configured model name, then that "
         "model will be used.  Otherwise, the given value will be treated as a file name and Longbow will attempt to "
         "read in the file and create a LibraryModel from it.  Longbow will assume the contents are the configuration "
         "of a LibraryModel as per LibraryModel.to_json()."
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
    "-l",
    "--min-length",
    type=int,
    default=0,
    show_default=True,
    required=False,
    help="Minimum length of a read to process.  Reads shorter than this length will not be annotated."
)
@click.option(
    "-L",
    "--max-length",
    type=int,
    default=longbow.utils.constants.DEFAULT_MAX_READ_LENGTH,
    show_default=True,
    required=False,
    help="Maximum length of a read to process.  Reads longer than this length will not be annotated."
)
@click.option(
    "--min-rq",
    type=float,
    default=-2,
    show_default=True,
    required=False,
    help="Minimum ccs-determined read quality for a read to be annotated.  CCS read quality range is [-1,1]."
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
def main(sirv_fasta, pbi, threads, output_file, model, chunk, min_length, max_length, min_rq, force, input_bam):
    """Prepare a dataset for use with model training."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    # Check to see if the output files exist:
    bam_utils.check_for_preexisting_files(output_file, exist_ok=force)

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

    # Load the SIRV reference sequences
    sirvs = load_fasta(sirv_fasta)

    # Get our model:
    if LibraryModel.has_prebuilt_model(model):
        m = LibraryModel.build_pre_configured_model(model)
    else:
        logger.info(f"Loading model from json file: %s", model)
        m = LibraryModel.from_json_file(model)
    logger.info(f"Using %s: %s", model, m.description)

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

            logger.info("Extracting training elements in %d reads from chunk %d/%d (reads %d-%d)", read_count, chunk, num_chunks, read_num, read_num + read_count - 1)
        else:
            read_count = bam_utils.load_read_count(pbi)
            logger.info("Extracting training elements in %d reads", read_count)

    # Create queues for data:
    queue_size = threads * 2 if threads < 10 else 20
    manager = mp.Manager()
    input_data_queue = manager.Queue(maxsize=queue_size)
    results = manager.Queue()

    # Start worker sub-processes:
    worker_pool = []

    for i in range(threads):
        p = mp.Process(
            target=_worker_segmentation_fn, args=(input_data_queue, results, i, m, min_length, max_length, min_rq, sirvs)
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

        # Get our header from the input bam file:
        out_header = bam_utils.create_bam_header_with_program_group(logger.name, bam_file.header, models=[m])

        # Start output worker:
        res = manager.dict({"num_reads_annotated": 0, "num_sections": 0})
        output_worker = mp.Process(
            target=_write_thread_fn,
            args=(results, output_file, not sys.stdin.isatty(), res, read_count, m)
        )
        output_worker.start()

        # Add in a sentinel value at the end of the queue - one for each subprocess - so we guarantee
        # that all subprocesses will exit:
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

    logger.info(
        f"Annotated {res['num_reads_annotated']} reads with {res['num_sections']} total sections."
    )
    et = time.time()
    logger.info(f"Done. Elapsed time: {et - t_start:2.2f}s. "
                f"Overall processing rate: {res['num_reads_annotated']/(et - t_start):2.2f} reads/s.")


def load_fasta(fasta):
    fa = FastaFile(fasta)
    return {key: value for (key, value) in list(map(lambda x: (x, fa.fetch(x)), fa.references))}


def get_segments(read):
    """Get the segments corresponding to a particular read by reading the segments tag information."""
    return read.to_string(), [
        SegmentInfo.from_tag(s) for s in read.get_tag(longbow.utils.constants.SEGMENTS_TAG).split(
            longbow.utils.constants.SEGMENT_TAG_DELIMITER)
    ]


def _write_thread_fn(out_queue, out_file_name, disable_pbar, res, read_count, model):
    """Thread / process fn to write out all our data."""

    with gzip.open(out_file_name, "wt") as out_file, tqdm.tqdm(
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
            rn = raw_data['rn']
            rq = raw_data['rq']
            example_sequences = raw_data['ex']

            # Write our training file:
            for name in example_sequences:
                rs = example_sequences[name]
                for r in rs:
                    out_file.write("\t".join([
                        rn,
                        name,
                        str(rq),
                        str(r.mismatch_count),
                        str(r.insertion_count),
                        str(r.deletion_count),
                        r.reference,
                        r.query,
                        r.alignment[0],
                        r.alignment[1],
                        r.alignment[2],
                        "\n"
                    ]))

            # Increment our counters:
            res["num_reads_annotated"] += 1 if len(example_sequences) > 0 else 0
            res["num_sections"] += len(example_sequences)

            pbar.update(1)


def _worker_segmentation_fn(in_queue, out_queue, worker_num, model, min_length, max_length, min_rq, sirvs):
    """Function to run in each subthread / subprocess.
    Segments each read and place the segments in the output queue."""

    num_reads_segmented = 0

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

        # Check for min/max length and min quality:

        if len(read.query_sequence) < min_length:
            logger.debug(f"Read is shorter than min length.  "
                         f"Skipping: {read.query_name} ({len(read.query_sequence)} < {min_length})")
            continue
        elif len(read.query_sequence) > max_length:
            logger.debug(f"Read is longer than max length.  "
                         f"Skipping: {read.query_name} ({len(read.query_sequence)} > {max_length})")
            continue
        elif read.get_tag("rq") < min_rq:
            logger.debug(f"Read quality is below the minimum.  "
                         f"Skipping: {read.query_name} ({read.get_tag('rq')} < {min_rq})")
            continue

        # Process and place our data on the output queue:
        segment_info = _segment_read(read, model)

        example_sequences = _extract_training_sequences(read, model, segment_info, sirvs)

        out_queue.put({
            'rn': read.query_name,
            'rq': read.get_tag('rq') if read.has_tag('rq') else -1,
            'ex': example_sequences
        })
        num_reads_segmented += 1

    logger.debug(f"Worker %d: Num reads segmented: %d", worker_num, num_reads_segmented)


def _extract_training_sequences(read, model, segment_info, sirvs):
    i = -1
    label = ""
    ppath = segment_info[1]
    query_sequence = read.query_sequence if segment_info[2] == False else bam_utils.reverse_complement(read.query_sequence)

    ssw_aligner = ssw.Aligner()
    example_sequences = {}

    while i < len(ppath) - 1:
        i += 1

        if ppath[i] != label:
            label = ppath[i]

            idx_left = i
            idx_right = i
            while idx_right < len(ppath) and ppath[idx_right] == label:
                idx_right += 1

            i = idx_right + 1

            q = query_sequence[idx_left:idx_right]
            if label in model.adapter_dict:
                if label not in ['CBC', 'UMI', 'Poly_A', 'VENUS', 'BOREAS', 'MARS']:
                    if label not in example_sequences:
                        example_sequences[label] = []

                    q = query_sequence[idx_left:idx_right]

                    refs = {}
                    if label == 'cDNA':
                        ref_names, _ = _prioritize_by_jaccard_index(q, sirvs, n=3)
                        for sirv_name in ref_names:
                            refs[sirv_name] = sirvs[sirv_name]

                    else:
                        refs[label] = model.adapter_dict[label]

                    best_score, best_ref_name, best_ref, best_r = 0, None, None, None
                    for ref_name in refs:
                        r = ssw_aligner.align(query=q, reference=refs[ref_name], revcomp=False)
                        if r.score > best_score:
                            best_score = r.score
                            best_ref_name = ref_name
                            best_ref = refs[ref_name]
                            best_r = r

                    if best_r.mismatch_count + best_r.deletion_count + best_r.insertion_count < int(0.15*len(q)):
                        r2, new_idx_left, new_idx_right = _scan_for_best_alignment(best_ref, query_sequence, idx_left, idx_right)

                        # print(f'{label} {refs[label] if label != "cDNA" else ""} {idx_left} {idx_right} {new_idx_left} {new_idx_right} {best_r.reference_coverage} {r2.reference_coverage} {best_r.query_coverage} {r2.query_coverage}')
                        # print(best_r.alignment_report())
                        # print(r2.alignment_report())
                        # print("")

                        if r2.reference_begin == 0 and r2.query_begin == 0 and r2.reference_coverage >= 0.95 and r2.query_coverage >= 0.95:
                            example_sequences[label].append(r2)


    return example_sequences


def _scan_for_best_alignment(ref_sequence, query_sequence, idx_left, idx_right):
    ssw_aligner = ssw.Aligner()

    best_coverage_sum_left, best_idx_left = 0, 0
    for shift_left in range(-min(200, int(len(ref_sequence)/5)), min(200, int(len(ref_sequence)/5))):
        if idx_left + shift_left >= 0 and idx_left + shift_left + 1 < len(query_sequence):
            q = query_sequence[max(0, idx_left + shift_left):max(idx_left + shift_left + 1, idx_right)]
            r = ssw_aligner.align(query=q, reference=ref_sequence, revcomp=False)

            if r.reference_coverage + r.query_coverage > best_coverage_sum_left:
                best_coverage_sum_left = r.reference_coverage + r.query_coverage
                best_idx_left = idx_left + shift_left

    best_coverage_sum_right, best_idx_right = 0, 0
    best_r = None
    for shift_right in range(-min(200, int(len(ref_sequence)/5)), min(200, int(len(ref_sequence)/5))):
        if idx_right + shift_right < len(query_sequence):
            q = query_sequence[best_idx_left:min(max(best_idx_left + 1, idx_right + shift_right), len(query_sequence))]
            r = ssw_aligner.align(query=q, reference=ref_sequence, revcomp=False)

            if r.reference_coverage + r.query_coverage > best_coverage_sum_right:
                best_coverage_sum_right = r.reference_coverage + r.query_coverage
                best_idx_right = idx_right + shift_right
                best_r = r

    return best_r, best_idx_left, best_idx_right


def _prioritize_by_jaccard_index(query, ref, k=7, n=3):
    jaccards = {}
    seqs = {}

    query_kmers = _kmerize(query, k)
    for name in ref:
        ref_kmers = _kmerize(ref[name], k)

        intersection = len(query_kmers.intersection(ref_kmers))
        union = len(query_kmers.union(ref_kmers))

        ji = intersection / union
        jaccards[name] = ji
        seqs[name] = ref[name]

    jaccards_sorted = {k: v for k, v in sorted(jaccards.items(), key=lambda item: item[1], reverse=True)}

    return list(jaccards_sorted.keys())[0:n], jaccards_sorted


def _kmerize(seq, k):
    kmers = set()

    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        kmers.add(kmer)

    return kmers


def _segment_read(read, model):
    is_rc = False
    logp, ppath = model.annotate(read.query_sequence)

    rc_logp, rc_ppath = model.annotate(bam_utils.reverse_complement(read.query_sequence))

    # print(f"Forward Path: {logp}: {ppath}")
    # print(f"Reverse Path: {rc_logp}: {rc_ppath}")

    if rc_logp > logp:
        logp = rc_logp
        ppath = rc_ppath
        is_rc = True
        logger.debug("Sequence scored better in RC: %s", read.query_name)

    return read.to_string(), ppath, logp, is_rc
