#!/usr/bin/env python

import os
import pathlib

import numpy as np
import pandas as pd
from plotnine import aes, geom_raster, ggplot, scale_fill_gradient


def draw_heatmap(CACHE_DIR: os.PathLike, target_series: str, ext: int) -> None:
    df_meta = pd.read_csv(CACHE_DIR / "meta.csv")

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
        return

    x_columns = [str(x) for x in range(-ext, ext + 1)]
    data = (
        pd.concat(ontargets)
        .groupby("sgRNA")[x_columns]
        .sum()
        .assign(total=lambda df: df.sum(axis=1))
        .sort_values("total")
        .reset_index()
        .reset_index(names="y")
        .melt(id_vars="y", value_vars=x_columns, var_name="x", value_name="score")
        .astype({"x": int})
    )

    color_max = np.percentile(data.loc[data["x"] == 0, "score"], q=95)

    (
        ggplot(data, mapping=aes(x="x", y="y", fill="score"))
        + geom_raster()
        + scale_fill_gradient(low="#FFFFFF", high="#FF0000", limits=[0, color_max])
    ).save(CACHE_DIR / "result" / f"{target_series}.pdf")


for target_series in [
    "hiplex1",
    "hiplex2",
    "hiplex3",
    "48sgrna",
]:
    draw_heatmap(
        CACHE_DIR=pathlib.Path(os.environ["CACHE_DIR"]),
        target_series=target_series,
        ext=17,
    )
