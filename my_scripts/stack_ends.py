#!/usr/bin/env python

import os
import pathlib

import pandas as pd
from plotnine import aes, geom_tile, ggplot, scale_fill_manual, scale_x_continuous


def stack_ends(target_series: str):
    df = (
        pd
        .read_csv(
            CACHE_DIR / "result" / "context_agg" / f"{target_series}.csv", header=0
        )
        .query("dominant_pos == -2")
        .reset_index(drop=True)
    )
    context = df.loc[0, "context"]
    df = df.query("context == @context and count > 0").reset_index(drop=True)[
        ["context", "pos", "count"]
    ]
    df = (
        pd
        .concat(
            [
                df,
                df["context"]
                .str.split("", expand=True)
                .drop(columns=[0, 35])
                .rename(columns={i: i - 18 for i in range(1, 35)}),
            ],
            axis=1,
        )
        .drop(columns=["context"])
        .melt(
            id_vars=["pos", "count"],
            value_vars=list(range(-17, 17)),
            var_name="x",
            value_name="nuc",
        )
        .astype({"x": int})
        .assign(x=lambda df: df["x"] + 0.5)
        .query("x > pos")
        .sort_values(by=["pos", "count"], ascending=False)
        .reset_index(drop=True)
        .assign(id=lambda df: pd.factorize(df["pos"])[0])
    )

    height = df.groupby("id")["count"].first()
    top = height.cumsum()
    bottom = [0] + top[:3].to_list()
    y = (top + bottom) / 2
    df = df.assign(
        width=1,
        height=lambda df: df["id"].map(height),
        y=lambda df: df["id"].map(y),
    )

    pdf_file = CACHE_DIR / "result" / "stack_ends" / f"{target_series}.pdf"
    pdf_file.parent.mkdir(parents=True, exist_ok=True)
    (
        ggplot(df, aes(x="x", y="y", fill="nuc"))
        + geom_tile(aes(width="width", height="height"))
        + scale_fill_manual(
            values={
                "A": "#0000FF",
                "C": "#FF0000",
                "G": "#00FF00",
                "T": "#FFFF00",
            }
        )
        + scale_x_continuous(limits=[-18, 18])
    ).save(pdf_file)


CACHE_DIR = pathlib.Path(os.environ["CACHE_DIR"])

stack_ends(target_series="hiplex1")
