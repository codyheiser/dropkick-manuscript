# -*- coding: utf-8 -*-
"""
Aggregate dropkick stats from multiple .h5ad files

@author: C Heiser
"""
import os, argparse
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from datetime import datetime
from distutils.util import strtobool
from sklearn.metrics import jaccard_score, roc_auc_score, roc_curve, auc
from progressbar import Percentage, ProgressBar, Bar, ETA


def set_diff(
    adata, labels=["dropkick_label", "CellRanger_label"], metrics=None, verbose=True
):
    """
    calculate sensitivity, specificity, and jaccard score between first and second labels

    Parameters:
        adata (anndata.AnnData): object with cell labels in .obs
        labels (list of str): two labels (columns of .obs) to compare the cell sets
            1 = real cell, 0 = empty or dead. labels[0] is prediction, labels[1] is ref.
        metrics (list of str): .obs metrics to average across barcodes for each label
        verbose (bool): print results to console or just return dict?

    Returns:
        out (dict): total and unique barcodes to each label, averages of each metric
            for each label
        prints results to console if verbose==True
    """
    if len(labels) != 2:
        raise ValueError("Please provide exactly two cell labels.")
    shared = len(
        set(
            adata.obs_names[adata.obs[labels[0]].isin(["True", True, 1.0, 1])]
        ).intersection(
            set(adata.obs_names[adata.obs[labels[1]].isin(["True", True, 1.0, 1])])
        )
    )
    unique_0 = len(
        set(
            adata.obs_names[adata.obs[labels[0]].isin(["True", True, 1.0, 1])]
        ).difference(
            set(adata.obs_names[adata.obs[labels[1]].isin(["True", True, 1.0, 1])])
        )
    )
    unique_1 = len(
        set(
            adata.obs_names[adata.obs[labels[1]].isin(["True", True, 1.0, 1])]
        ).difference(
            set(adata.obs_names[adata.obs[labels[0]].isin(["True", True, 1.0, 1])])
        )
    )
    # initiate dict for output
    out = {}

    # add total and unique barcodes for labels[0]
    out["{}_total".format(labels[0])] = (
        adata.obs[labels[0]].isin(["True", True, 1.0, 1]).sum()
    )
    # print results for labels[0] to console
    if verbose:
        print(
            "{} cells in {} - {} unique".format(
                adata.obs[labels[0]].isin(["True", True, 1.0, 1]).sum(),
                labels[0],
                unique_0,
            )
        )
    if metrics is not None:
        if isinstance(metrics, str):
            metrics = [metrics]
        for m in metrics:
            out["{}_{}_avg".format(labels[0], m)] = adata.obs.loc[
                adata.obs[labels[0]].isin(["True", True, 1.0, 1]), m
            ].mean()
            if verbose:
                print(
                    "\t{}: {:0.3}".format(
                        m,
                        round(
                            adata.obs.loc[
                                adata.obs[labels[0]].isin(["True", True, 1.0, 1]), m
                            ].mean(),
                            3,
                        ),
                    )
                )

    # add total and unique barcodes for labels[1]
    out["{}_total".format(labels[1])] = (
        adata.obs[labels[1]].isin(["True", True, 1.0, 1]).sum()
    )
    # print results for labels[1] to console
    if verbose:
        print(
            "{} cells in {} - {} unique".format(
                adata.obs[labels[1]].isin(["True", True, 1.0, 1]).sum(),
                labels[1],
                unique_1,
            )
        )
    if metrics is not None:
        for m in metrics:
            out["{}_{}_avg".format(labels[1], m)] = adata.obs.loc[
                adata.obs[labels[1]].isin(["True", True, 1.0, 1]), m
            ].mean()
            if verbose:
                print(
                    "\t{}: {:0.3}".format(
                        m,
                        round(
                            adata.obs.loc[
                                adata.obs[labels[1]].isin(["True", True, 1.0, 1]), m
                            ].mean(),
                            3,
                        ),
                    )
                )

    out["{}_sensitivity".format(labels[1])] = (
        shared / adata.obs[labels[1]].isin(["True", True, 1.0, 1]).sum()
    )
    out["{}_specificity".format(labels[1])] = (
        adata.n_obs - adata.obs[labels[1]].isin(["True", True, 1.0, 1]).sum()
    ) / (
        adata.n_obs - adata.obs[labels[1]].isin(["True", True, 1.0, 1]).sum() + unique_0
    )

    # ensure workable dtype for jaccard scoring
    if "True" in adata.obs[labels[0]].unique():
        adata.obs[labels[0]] = adata.obs[labels[0]].apply(strtobool)
    if "True" in adata.obs[labels[1]].unique():
        adata.obs[labels[1]] = adata.obs[labels[1]].apply(strtobool)

    out["{}_jaccard".format(labels[1])] = jaccard_score(
        adata.obs[labels[0]], adata.obs[labels[1]]
    )

    return out


def agg_stats(adatas, labels, names, compare_label="dropkick_label", compare_score="dropkick_score", roc=True, outdir="."):
    out = None
    for l in labels:
        if roc:
            plt.figure(figsize=(7, 5))
            plt.plot([0, 1], [0, 1], color="k", lw=2, linestyle="--")
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel("False Positive Rate", fontsize=12)
            plt.ylabel("True Positive Rate", fontsize=12)
            plt.tick_params(labelsize=12)
        df = pd.DataFrame()
        tmp = {}
        widgets = [
            'Calculating stats against "{}": '.format(l),
            Percentage(),
            " ",
            Bar(marker="#", left="[", right="]"),
            " ",
            ETA(),
        ]
        pbar = ProgressBar(widgets=widgets, maxval=len(adatas))
        pbar.start()
        i = 0
        for a, n in zip(adatas, names):
            tmp = {**tmp, **set_diff(a, labels=[compare_label, l], verbose=False)}
            if roc:
                fpr, tpr, _ = roc_curve(a.obs[l], a.obs[compare_score])
                tmp["{}_auc".format(l)] = auc(fpr, tpr)
                plt.plot(
                    fpr,
                    tpr,
                    lw=1.5,
                    alpha=0.7,
                    label=n, #"{} (AUC = {})".format(n, round(tmp["{}_auc".format(l)], 3)),
                )
            df = pd.concat([df, pd.DataFrame(tmp, index=[n])])
            i += 1
            pbar.update(i)
        pbar.finish()

        if out is None:
            out = df.copy()
        else:
            out = out.merge(
                df,
                how="outer",
                left_index=True,
                right_index=True,
                suffixes=(None, "_x"),
            )

        if roc:
            plt.legend(fontsize=10, bbox_to_anchor=(1.05, 1.0), loc="upper left")
            plt.tight_layout()
            print('Saving ROC plot for "{}" to roc_{}.png'.format(l, l))
            plt.savefig("{}/roc_{}.png".format(outdir, l))

    if "{}_total_x".format(compare_label) in out.columns:
        out.drop(columns=["{}_total_x".format(compare_label)], inplace=True)
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "files",
        type=str,
        nargs="*",
        help="List of paths to .h5ad files to aggregate stats for.",
    )
    parser.add_argument(
        "-l",
        "--labels",
        type=str,
        help="Labels to compare compare_label to (i.e. CellRanger_3)",
        nargs="*",
        default=["CellRanger_2", "CellRanger_3"],
    )
    parser.add_argument(
        "-cl",
        "--compare-label",
        type=str,
        help="Label to compare (default dropkick_label)",
        default="dropkick_label",
    )
    parser.add_argument(
        "-cs",
        "--compare-score",
        type=str,
        help="Score to compare (default dropkick_score)",
        default="dropkick_score",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        help="Path to output directory. Default './'",
        nargs="?",
        default="./",
    )
    args = parser.parse_args()

    # read the .h5ad files into list
    adatas = []
    print("\nReading: ", end="")
    for f in args.files:
        # read file into anndata obj
        print("{} ".format(f), end="")
        adatas.append(sc.read(f))
    print("\n")

    # aggregate stats into pd.DataFrame
    df = agg_stats(
        adatas,
        labels=args.labels,
        names=[os.path.splitext(os.path.basename(x))[0] for x in args.files],
        compare_label=args.compare_label,
        compare_score=args.compare_score,
        roc=True,
        outdir=args.outdir,
    )

    # save output df to .csv file
    print("Writing stats to {}/stats_{}.csv".format(args.outdir, "_".join(args.labels)))
    df.to_csv("{}/stats_{}.csv".format(args.outdir, "_".join(args.labels)))
