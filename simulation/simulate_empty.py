import os, errno, argparse
import numpy as np
import scanpy as sc
from numpy.random import default_rng


def check_dir_exists(path):
    """
    Checks if directory already exists or not and creates it if it doesn't

    Parameters:
        path (str): path to directory

    Returns:
        tries to make directory at path, unless it already exists
    """
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def partition_adata(
    adata,
    n_top_cells=None,
    n_bottom_counts=100,
    n_ambient=1000,
    min_genes=5,
    min_cells=5,
):
    """
    Splits anndata into ground-truth cells and empty droplets by barcode ranking,
        performs pre-filtering of cells and genes, and defines top ambient genes
    """
    adata.obs["n_counts"] = adata.X.sum(axis=1)
    adata.obs["partition"] = ""
    adata.obs.loc[adata.obs.n_counts <= n_bottom_counts, "partition"] = "bottom"
    if isinstance(n_top_cells, int):
        adata.obs.loc[
            adata.obs.n_counts.nlargest(n_top_cells, keep="all").index, "partition"
        ] = "top"
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata.var["pct_dropout_by_counts"] = np.array(
        (1 - (adata.X.astype(bool).sum(axis=0) / adata.n_obs)) * 100
    ).squeeze()
    adata.var["top_ambient"] = "False"
    adata.var.loc[
        adata.var.pct_dropout_by_counts.nsmallest(n_ambient).index, "top_ambient"
    ] = "True"
    return adata


def simulate_empties(adata, seed=18, n_empties=2000, min_counts=10, max_counts=100):
    """
    Simulates empty droplets based on partitions from partition_adata(),
        using random multinomial number generator
    """
    mtx = (
        adata[adata.obs.partition == "bottom", adata.var.top_ambient == "True"]
        .X.copy()
        .todense()
    )
    sums = np.array(mtx.sum(axis=0)).squeeze()
    sums = sums / sums.sum()  # get probs that add up to 1.0
    rng = default_rng(seed)  # initiate number generator
    empties = rng.multinomial(
        n=rng.integers(low=min_counts, high=max_counts, size=n_empties),
        pvals=sums,
        size=n_empties,
    )
    adata_e = sc.AnnData(
        empties, var=adata.var.loc[adata.var.top_ambient == "True", []]
    )
    return adata_e


def run_sim(args):
    """
    Builds simulated dataset from ground-truth real cells and simulated
        empty droplets, writes results to .h5ad file
    """
    # get basename of file for writing outputs
    name = os.path.splitext(os.path.basename(args.file))[0]
    # read file into anndata obj
    print("Reading {}".format(args.file), end="")
    a = sc.read(args.file)
    # print information about counts
    print(" - {} cells and {} genes".format(a.shape[0], a.shape[1]))
    # partition data
    print("Partitioning counts data...")
    a = partition_adata(
        a,
        n_top_cells=args.n_top_cells,
        n_bottom_counts=args.n_bottom_counts,
        n_ambient=args.n_ambient,
    )
    # simulate empty droplets
    print("Simulating empty droplets...")
    a_e = simulate_empties(
        a,
        seed=args.seed,
        n_empties=args.n_empties,
        min_counts=args.min_counts_empty,
        max_counts=args.max_counts_empty,
    )
    if isinstance(args.n_top_cells, int):
        # add simulated empties to ground-truth cells only
        a_sim = a[a.obs.partition == "top", :].copy()
        a_sim = a_sim.concatenate(
            a_e,
            join="outer",
            batch_key="simulated",
            batch_categories=["False", "True"],
            fill_value=0,
        )
    else:
        # add simulated empties to whole partitioned counts matrix
        a_sim = a.concatenate(
            a_e,
            join="outer",
            batch_key="simulated",
            batch_categories=["False","True"],
            fill_value=0,
        )
    # save file as .h5ad
    print("Writing simulation to {}/{}.h5ad".format(args.outdir, name))
    check_dir_exists(args.outdir)
    a_sim.write("{}/{}_sim.h5ad".format(args.outdir, name), compression="gzip")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "file",
        type=str,
        help="Path to .h5ad file",
    )
    parser.add_argument(
        "--n-top-cells",
        type=int,
        help="Number of top barcodes by counts to use as ground-truth cells. If None, keep all barcodes from original matrix. Default None.",
        default=None,
    )
    parser.add_argument(
        "--n-bottom-counts",
        type=int,
        help="Total counts below which to label barcodes as ground-truth empty droplets. Default 100.",
        default=100,
    )
    parser.add_argument(
        "--n-ambient",
        type=int,
        help="Number of bottom genes by dropout rate to use as ground-truth ambient genes. Default 1000.",
        default=1000,
    )
    parser.add_argument(
        "--n-empties",
        type=int,
        help="Number of empty droplets to simulate. Default 2000.",
        default=2000,
    )
    parser.add_argument(
        "--min-counts-empty",
        type=int,
        help="Minimum total counts in simulated empties. Default 10.",
        default=10,
    )
    parser.add_argument(
        "--max-counts-empty",
        type=int,
        help="Maximum total counts in simulated empties. Default 100.",
        default=100,
    )
    parser.add_argument(
        "-s",
        "--seed",
        type=int,
        help="Random number generator seed. Default 18.",
        default=18,
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        nargs="?",
        default=".",
        help="Output directory for writing h5ad. Default './'",
    )

    parser.set_defaults(func=run_sim)
    args = parser.parse_args()
    args.func(args)
