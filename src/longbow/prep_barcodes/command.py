import logging
import sys
import time

import click
import click_log

from ..utils import barcode_utils


logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("prep_barcodes")
click_log.basic_config(logger)


@click.command(name=logger.name)
@click_log.simple_verbosity_option(logger)
@click.option(
    "-o",
    "--output-barcodes",
    type=str,
    help="name to give to output file",
)
@click.option(
    "-k",
    "--kmer-size",
    default=3,
    show_default=True,
    help="Kmer size to use when creating barcode vectors."
)
@click.option(
    "-n",
    "--n-neighbors",
    default=10,
    show_default=True,
    help="Number of neighbors to use by default for kneighbors queries."
)
#@click.argument("input-barcodes", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
@click.argument("input-barcodes", type=str)
def main(output_barcodes, kmer_size, n_neighbors, input_barcodes):
    """Prepare barcode database for use with subsequent error correction."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    logger.info(f"Reading barcodes...")
    barcodes = barcode_utils.load_barcode_whitelist(input_barcodes)

    logger.info(f"Creating MIPS database (k={kmer_size}, n_neighbors={n_neighbors})...")
    knn = barcode_utils.prepare_database(barcodes, k=kmer_size, n_neighbors=n_neighbors)

    logger.info(f"Saving database to {output_barcodes}...")
    barcode_utils.save_database(knn, output_barcodes)

    et = time.time()
    logger.info(f"Done. Elapsed time: {et - t_start:2.2f}s.")