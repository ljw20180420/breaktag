#!/usr/bin/env python

import os
import pathlib
import sys
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Seq
from plotnine import aes, geom_raster, ggplot, scale_fill_gradient2
from roukos import load_fasta


def get_context(row: pd.Series, genome: OrderedDict) -> str:
    chromsome = genome[row["chrom"]]
    context = chromsome.ff.fetch(
        chromsome.ref, row["cut"] - 17, row["cut"] + 17
    ).upper()
    if row["sense"] == "-":
        context = str(Seq.Seq(context).reverse_complement())
    return context


def collect_data(target_series: str) -> pd.DataFrame:
    df_meta = pd.read_csv(CACHE_DIR / "meta.csv", header=0)

    genome = load_fasta(
        CACHE_DIR
        / "human_dataset"
        / "ncbi_dataset"
        / "data"
        / "GCF_000001405.40"
        / "GCF_000001405.40_GRCh38.p14_genomic.fna"
    )

    ontargets = []
    for sample_title, series, target, nontarget, sgRNA_lib in zip(
        df_meta["Sample_title"],
        df_meta["series"],
        df_meta["target"],
        df_meta["nontarget"],
        df_meta["sgRNA_lib"],
    ):
        if nontarget == "unknown" or sgRNA_lib == "unknown" or series != target_series:
            continue

        ontargets.append(
            pd.read_csv(CACHE_DIR / "result" / f"{sample_title}.{series}.csv").query(
                'sgRNA == sgRNA_lib and tnt == "target" and PAM.str.endswith("GG")'
            )
        )

    if len(ontargets) == 0:
        print(f"{target_series} has empty results", file=sys.stderr)
        return

    return (
        pd
        .concat(ontargets)
        .reset_index(drop=True)
        .assign(
            context=lambda df: df.apply(
                lambda row, genome=genome: get_context(row, genome), axis=1
            )
        )
    )


def get_dominant_pos(df: pd.DataFrame) -> pd.Series:
    df["pos"] = df.loc[df["count"].idxmax(), "pos"]
    return df["pos"]


def context_agg(df: pd.DataFrame, ext: int) -> pd.DataFrame:
    x_columns = [str(x) for x in range(-ext, ext + 1)]
    df = (
        df
        .assign(down=lambda df: df["strand"] == df["sense"])
        .groupby(["context", "down"])[x_columns]
        .sum()
        .reset_index()
        .melt(
            id_vars=["context", "down"],
            value_vars=x_columns,
            var_name="pos",
            value_name="count",
        )
        .astype({"pos": int})
        .assign(
            dominant_pos=lambda df: df.groupby(
                ["context", "down"], group_keys=False
            ).apply(get_dominant_pos, include_groups=False),
            dominant_count=lambda df: df.groupby(["context", "down"], group_keys=False)[
                "count"
            ].transform("max"),
            total_count=lambda df: df.groupby(["context", "down"], group_keys=False)[
                "count"
            ].transform("sum"),
            dominant_ratio=lambda df: df["dominant_count"] / df["total_count"],
            ratio=lambda df: df["count"] / df["total_count"],
        )
    )

    df_down = df.query("down").reset_index(drop=True)
    df_up = df.query("not down").reset_index(drop=True)
    for direc, df_direc in zip(["down", "up"], [df_down, df_up]):
        df_direc = (
            df_direc
            .sort_values(
                by=[
                    "dominant_pos",
                    "dominant_ratio",
                    "context",
                ]
            )
            .reset_index(drop=True)
            .assign(id=lambda df: pd.factorize(df["context"])[0])
        )

        csv_file = CACHE_DIR / "result" / "context_agg" / f"{target_series}_{direc}.csv"
        csv_file.parent.mkdir(parents=True, exist_ok=True)
        df_direc.to_csv(csv_file, index=False)

        pdf_file = CACHE_DIR / "result" / "context_agg" / f"{target_series}_{direc}.pdf"
        pdf_file.parent.mkdir(parents=True, exist_ok=True)
        color_max = np.percentile(df_direc.loc[df_direc["pos"] == 0, "ratio"], q=100)
        (
            ggplot(df_direc, mapping=aes(x="pos", y="id", fill="ratio"))
            + geom_raster(interpolate="nearest")
            + scale_fill_gradient2(
                low="#0707D6",
                mid="#FFFFFF",
                high="#FF0000",
                midpoint=color_max / 2,
                limits=[0, color_max],
            )
        ).save(pdf_file)
        plt.close("all")


if __name__ == "__main__":
    CACHE_DIR = pathlib.Path(os.environ["CACHE_DIR"])
    target_series = "hiplex1"
    df = collect_data(target_series)
    context_agg(df, ext=17)
