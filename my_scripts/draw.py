#!/usr/bin/env python

import os
import pathlib

import pandas as pd
from plotnine import aes, geom_raster, ggplot

CACHE_DIR = pathlib.Path(os.environ["CACHE_DIR"])
ext = 17

df_meta = pd.read_csv(CACHE_DIR / "meta.csv")

ontargets = []
for sample_title, series, target, nontarget, sgRNA_lib in zip(
    df_meta["Sample_title"],
    df_meta["series"],
    df_meta["target"],
    df_meta["nontarget"],
    df_meta["sgRNA_lib"],
):
    if nontarget == "unknown" or sgRNA_lib == "unknown" or series != "hiplex1":
        continue

    ontargets.append(
        pd.read_csv(CACHE_DIR / "result" / f"{sample_title}.{series}.csv").query(
            'sgRNA == sgRNA_lib and tnt == "target" and PAM.str.endswith("GG")'
        )
    )

ontarget = pd.concat(ontargets).reset_index(drop=True)
x_columns = [str(x) for x in range(-ext, ext + 1)]
data = (
    ontarget[x_columns]
    .reset_index(names="y")
    .melt(id_vars="y", value_vars=x_columns, var_name="x", value_name="score")
    .astype({"x": int})
)


(ggplot(data, mapping=aes(x="x", y="y", fill="score")) + geom_raster()).save(
    CACHE_DIR / "result" / f"{sample_title}.{series}.pdf"
)
