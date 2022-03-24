import os.path
from pathlib import Path

CACHED_GENES_FILENAME = "cached_genes.tsv"
TENSOR_PATTERN = "tensors/tensor_gene_{}.npy"
PRED_CACHE_PATTERN = "pred_{}_{}_{}.npz"


def get_cached_genes_filename(data_folder):
    return os.path.join(data_folder, CACHED_GENES_FILENAME)


def get_tensor_filename(gene, data_folder):
    return os.path.join(data_folder, TENSOR_PATTERN.format(gene.id))


def get_pred_cache_filename(pred_cache_folder, chrom, start, stop):
    folder = os.path.join(pred_cache_folder, f"{chrom}/{str(start)[:3]}")
    Path(folder).mkdir(parents=True, exist_ok=True)
    return os.path.join(folder, PRED_CACHE_PATTERN.format(chrom, start, stop))
