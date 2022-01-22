import logging
import time
import math
import sys
import gzip
import re
import random
import pprint

import click
import click_log
import tqdm

import multiprocessing as mp
import concurrent.futures

import pysam
import pandas as pd
import numpy as np

import longbow.utils.constants
from ..utils import bam_utils, barcode_utils
from ..utils import model as LongbowModel
from ..utils.model import LibraryModel

from ..utils import bam_utils

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("train")
click_log.basic_config(logger)


@click.command(name="train")
@click_log.simple_verbosity_option(logger)
# @click.option(
#     "-s",
#     "--sirv-fasta",
#     required=True,
#     type=click.Path(),
#     help="SIRV fasta file",
# )
# @click.option(
#     "-n",
#     "--num-training-samples",
#     type=int,
#     default=10,
#     show_default=True,
#     help="number of training samples to use",
# )
# @click.option(
#     "-i",
#     "--max-training-iterations",
#     type=int,
#     default=5,
#     show_default=True,
#     help="number of training iterations to use",
# )
@click.option(
    "-t",
    "--threads",
    type=int,
    default=mp.cpu_count() - 1,
    show_default=True,
    help="number of threads to use (0 for all)",
)
# @click.option(
#     "-o",
#     "--output-yaml",
#     required=True,
#     type=click.Path(exists=False),
#     help="trained model",
# )
@click.option(
    "-m",
    "--model",
    default=longbow.utils.constants.DEFAULT_MODEL,
    show_default=True,
    help="The model to train.  If the given value is a pre-configured model name, then that "
         "model will be used.  Otherwise, the given value will be treated as a file name and Longbow will attempt to "
         "read in the file and create a LibraryModel from it.  Longbow will assume the contents are the configuration "
         "of a LibraryModel as per LibraryModel.to_json()."
)
@click.option(
    "-w",
    "--cbc-whitelist-db",
    required=False,
    type=click.Path(),
    help="Cell barcode whitelist.",
)
@click.argument("training-file", type=click.Path(exists=True))
#def main(sirv_fasta, num_training_samples, max_training_iterations, threads, output_yaml, model, training_bam):
def main(threads, model, cbc_whitelist_db, training_file):
    """Train transition and emission probabilities on real data."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

    # Load CBC whitelist database (empty set if no whitelist is specified)
    cbc_set = barcode_utils.load_barcode_whitelist(re.sub(".db$", "", cbc_whitelist_db))

    # Get our model:
    if LibraryModel.has_prebuilt_model(model):
        m = LibraryModel.build_pre_configured_model(model)
    else:
        logger.info(f"Loading model from json file: %s", model)
        m = LibraryModel.from_json_file(model)
    logger.info(f"Using %s: %s", model, m.description)

    # Load training data
    training_data = load_training_data(training_file)
    error_rate_tbl = compute_error_rates(training_data)

    # print(error_rate_tbl)

    simulate_reads(m, training_data, error_rate_tbl, cbc_set)

    # pp = pprint.PrettyPrinter(indent=4)
    # pp.pprint(training_data)

    # logger.info("Loaded %d training sequences", len(training_seqs))

    # logger.info("Starting training...", len(training_seqs))
    # improvement, history = m.fit(
    #     sequences=training_seqs,
    #     max_iterations=max_training_iterations,
    #     stop_threshold=1e-1,
    #     return_history=True,
    #     verbose=True,
    #     n_jobs=threads,
    # )

    # with open(output_yaml, "w") as model_file:
    #     print(improvement.to_yaml(), file=model_file)

    logger.info(f"Done. Elapsed time: %2.2fs.", time.time() - t_start)


def load_training_data(training_file):
    training_data = {}

    header = [ 'read_name', 'element', 'rq', 'mismatch_count', 'insertion_count', 'deletion_count', 'reference', 'query', 'align_reference', 'align_track', 'align_query' ]
    with gzip.open(training_file, 'rt') as tf:
        for l in tf:
            values = l.rstrip().split('\t')
            res = {header[i]: values[i] for i in range(len(header))}
            res['rq'] = float(res['rq'])
            res['mismatch_count'] = int(res['mismatch_count'])
            res['insertion_count'] = int(res['insertion_count'])
            res['deletion_count'] = int(res['deletion_count'])

            if res['element'] not in training_data:
                training_data[res['element']] = []

            training_data[res['element']].append(res)

    return training_data


def quantize_rq(rq, rq_bins = [-1.0, 0.0]):
    for i in range(len(rq_bins)):
        bin_left = rq_bins[i]
        bin_right = rq_bins[i+1] if i+1 < len(rq_bins) else 1.0

        if bin_left <= rq and rq < bin_right:
            return bin_left

    return -1.0


def compute_error_rates(training_data, rq_bins = [-1.0, 0.0]):
    tbl_header = ["rq", "mrate", "irate", "drate", "ins_len", "del_len"]
    tbl_rows = []

    for i in range(len(rq_bins)):
        bin_left = rq_bins[i]
        bin_right = rq_bins[i+1] if i+1 < len(rq_bins) else 1.0

        for e in training_data:
            for t in training_data[e]:
                if bin_left <= t['rq'] and t['rq'] < bin_right:
                    mrate = ( t['mismatch_count'] / len(t['reference']) ) + 0.0001
                    irate = ( t['insertion_count'] / len(t['reference']) ) + 0.0001
                    drate = ( t['deletion_count'] / len(t['reference']) ) + 0.0001
                    ins_len = None
                    del_len = None

                    if t['insertion_count'] > 0:
                        ins_len = np.mean(list(map(lambda x: len(x), re.findall("-+", t['align_reference']))))

                    if t['deletion_count'] > 0:
                        del_len = np.mean(list(map(lambda x: len(x), re.findall("-+", t['align_query']))))

                    tbl_rows.append([bin_left, mrate, irate, drate, ins_len, del_len])

    tbl_new = pd.DataFrame(tbl_rows, columns=tbl_header)

    return tbl_new.groupby('rq').mean()


def simulate_reads(m, training_data, error_rate_tbl, cbc_set):
    num_training_samples = {}
    training_data_iters = {}

    for rq in list(error_rate_tbl.index):
        num_training_samples[rq] = None
        training_data_iters[rq] = {}

    for k in training_data:
        for i in range(len(training_data[k])):
            rq = quantize_rq(training_data[k][i]['rq'])

            num_training_samples[rq] = min(list(filter(lambda x: x is not None, [len(training_data[k]), num_training_samples[rq]])))
            training_data_iters[rq][k] = iter(training_data[k])

    for rq in num_training_samples:
        n = num_training_samples[rq]

        for i in range(n):
            read = []

            idx1 = random.choice(list(range(len(m.array_element_structure))))
            idx2 = random.choice(list(range(len(m.array_element_structure))))
            while idx2 == idx1:
                idx2 = random.choice(list(range(len(m.array_element_structure))))

            for array_structures in m.array_element_structure[min(idx1, idx2):max(idx1, idx2)]:
                for k in array_structures:
                    seg = ''
                    if k not in training_data:
                        ref = m.adapter_dict[k]
                        if not isinstance(ref, dict):
                            seg = add_noise(ref, error_rate_tbl.loc[[rq]])
                        else:
                            if k == 'Poly_A':
                                rf1 = 'A'*int(random.uniform(27, 33))
                            elif k == 'UMI':
                                rf1 = simulate_random_sequence(ref['FixedLengthRandomBases'])
                            elif k == 'CBC':
                                rf1 = random.choice(cbc_set)

                            seg = add_noise(rf1, error_rate_tbl.loc[[rq]])
                    else:
                        seg = next(training_data_iters[rq][k])["query"]

                    read.append(seg)

            print(read)



def simulate_random_sequence(length):
    alphabet = ['A', 'C', 'G', 'T']
    
    bases = random.choices(alphabet, k=length)
    bases = ''.join(bases)
    
    return bases


def add_noise(seq, error_rate_tbl_rq):
    mismatch_rate = error_rate_tbl_rq['mrate'][0]
    insertion_rate = error_rate_tbl_rq['irate'][0]
    deletion_rate = error_rate_tbl_rq['drate'][0]
    ins_len = error_rate_tbl_rq['ins_len'][0]
    del_len = error_rate_tbl_rq['del_len'][0]

    mm_positions = int(mismatch_rate*len(seq))
    ins_positions = int(insertion_rate*len(seq))
    del_positions = int(deletion_rate*len(seq))
    
    newseq = list(seq)
    
    for i in range(mm_positions):
        p = int(random.uniform(0, len(seq)))
        newseq[p] = bam_utils.reverse_complement(newseq[p])

    for i in range(ins_positions):
        p = int(random.uniform(0, len(seq)))
        #l = int(random.uniform(2, 5))
        l = int(max(1.0, random.normal(loc=ins_len, scale=ins_len, size=1)))
        newseq[p] = l*newseq[p]
        
    for i in range(del_positions):
        p = int(random.uniform(0, len(seq)))
        newseq[p] = ''

    return ''.join(newseq)
    

def select_read(read, model):
    flogp = -math.inf
    fseq = None

    # Use the untrained model to determine if we should add this training
    # example in the forward or reverse-complement orientation.
    for seq in [read.query_sequence, bam_utils.reverse_complement(read.query_sequence)]:
        logp, ppath = model.annotate(seq, smooth_islands=True)

        if logp > flogp:
            flogp = logp
            fseq = seq

    return flogp, fseq
