import pandas as pd
import argparse
import logging
import os.path
from pathlib import Path

from args import add_args, setup
from filenames import get_cached_genes_filename


def load_genes(genes_filename, gene_type, data_folder):
    folder = os.path.join(data_folder, "tensors")
    Path(folder).mkdir(parents=True, exist_ok=True)

    cached_genes_filename = get_cached_genes_filename(data_folder)
    try:
        logging.info(f"Loading genes from {cached_genes_filename}...")
        genes = pd.read_csv(cached_genes_filename, sep="\t", index_col=0)
        logging.info(f"Loaded {len(genes):,} genes.")
    except FileNotFoundError:
        logging.info(f"Failed. Loading genes from: {genes_filename}...")
        genes = pd.read_csv(genes_filename)
        if gene_type == "pc":
            genes = genes[genes.type == "protein_coding"]
        elif gene_type == 'lincRNA':
            genes = genes[genes.type == "lincRNA"]
        elif gene_type == 'all':
            genes = genes[genes.type != "rRNA"]
        else:
            raise ValueError('gene_type has to be one of all, pc, and lincRNA')
        genes = genes[~genes.seqnames.isin(["chrX", "chrY"])]
        genes = genes.dropna()
        genes.TSS = genes.TSS.astype("int")
        genes.CAGE_representative_TSS = genes.CAGE_representative_TSS.astype("int")

        logging.info(f"Loaded {len(genes):,} genes with type={gene_type}.")

        genes.to_csv(cached_genes_filename, sep="\t")

    return genes


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    add_args(parser)
    args = parser.parse_args()
    setup(args)

    # load protein coding genes and Chromatin Contact Domains (CCDs), find CCDs containing genes
    genes, ccds = load_genes(args.genes_filename, args.gene_type, args.ccds_filename, args.data_folder)
