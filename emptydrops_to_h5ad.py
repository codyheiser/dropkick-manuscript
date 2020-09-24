import argparse
import pandas as pd
import scanpy as sc

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "adata",
        type=str,
        help="AnnData object saved as .h5ad file to add EmptyDrops labels to",
    )
    parser.add_argument(
        "emptydrops_out",
        type=str,
        help="Output from EmptyDrops as .csv file corresponding to adata",
    )
    args = parser.parse_args()

    a = sc.read(args.adata)  # read in anndata
    e = pd.read_csv(args.emptydrops_out, index_col=0)  # read in EmptyDrops results

    # add emptydrops results to anndata object in .obs
    a.obs[["EmptyDrops_LogProb","EmptyDrops_pval","EmptyDrops_limited","EmptyDrops_FDR"]] = e[["LogProb","PValue","Limited","FDR"]].values
    a.obs["EmptyDrops"] = "False"
    a.obs.loc[a.obs.EmptyDrops_FDR <= 0.001, "EmptyDrops"] = "True"

    # overwrite original .h5ad file
    a.write(args.adata, compression="gzip")
