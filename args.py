import logging
import numpy as np
import os
import platform

is_mbp = platform.node().startswith("MBP") or platform.node().startswith("MacBook")


def setup(args):
    # log all the SLURM environment variables if run on the SLURM node
    nodename = os.environ.get("SLURM_JOB_NODELIST", "")
    nodename = nodename + " " if nodename else ""
    logging.basicConfig(level=logging.INFO, format=nodename+'%(asctime)s %(message)s', datefmt='%H:%M:%S')

    logging.info("SLURM env:")
    for name in sorted(os.environ):
        if name.startswith("SLURM_"):
            logging.info(f"{name}: {os.environ[name]}")

    logging.info("Command line args:\n"+"\n".join([f"{k}: {v}" for k, v in args.__dict__.items()]))

    np.random.seed(1410)


def add_args(parser):
    # data folder
    parser.add_argument("--data_folder", type=str, default="data_beluga_expecto_hg38", help="Folder to store all the data")

    # genes
    parser.add_argument("--genes_filename", type=str,
                        default="./expecto/geneanno_GRCh38.csv", help="Genes csv file")
    parser.add_argument("--gene_type", type=str,
                        default="all", help="Gene type: pc, lincRNA, all (default:all)")

    # subset data for calculation parallelization
    parser.add_argument("--divisor", type=int, help="Divisor to parallelize computation")
    parser.add_argument("--reminder", type=int, help="Reminder to parallelize computation")

    # tensors
    parser.add_argument("--bin_size", type=int, default=200, help="Bin size (default: 200)")
    parser.add_argument("--window_size", type=int, default=2000, help="Window size (default: 2000)")
    parser.add_argument("--tensor_width", type=int, default=10, help="Tensor width (default: 10)")
    parser.add_argument("--features_count", type=int, default=2002, help="Features count (default: 2002)")
    parser.add_argument("--batch_size", type=int, default=32, help="Batch size (default: 32)")
    parser.add_argument("--genome_filename", type=str, default="./data/hg38.fa", help="Genome filename")
    parser.add_argument("--beluga_model_filename", type=str, default="./expecto/deepsea.beluga.pth",
                        help="Expecto Beluga model filename")
    parser.add_argument("--test_only", type=bool, default=False, help="Generate test tensors only (default: False)")

    # 3D args
    parser.add_argument("--mcool_filename", type=str, default="", help="HiC/CHiA PET file")
    parser.add_argument("--hic_resolution", type=int, default=1000, help="HiC resolution (default: 1000)")
    parser.add_argument('--max_dist', type=float, default=1.15, help="Maximum distance to calculate tensor")

    # Beluga prediction cache
    parser.add_argument("--pred_cache_folder", type=str, default="./pred_cache", help="Prediction cache folder")

    # CUDA
    parser.add_argument('--compute_platform', choices=["CPU", "CUDA"], default="CUDA", help="Computation platform (CPU or CUDA")

    # training
    parser.add_argument("--gene_expression_filename", type=str, default="./expecto/geneanno.exp.csv",
                        help="Gene expression filename")
    parser.add_argument('--pseudocount', type=float, default=0.0001, help="pseudocount")
    parser.add_argument('--target_index', type=int, default=22, help="Tissue Index (default: 22)")
    parser.add_argument("--biosample", type=str, default="", help="Biosample name")

    # logging
    parser.add_argument("--google_sheets_credentials_filename", type=str, default="",
                        help="Google Sheet credentials json file")
    parser.add_argument("--google_spreadsheet_id", type=str, default="",
                        help="Google Sheet Id")
    parser.add_argument("--log_note", type=str, default="Beluga Expecto Hg38+HiC GM12878", help="Log Note")
