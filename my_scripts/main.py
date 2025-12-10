#!/usr/bin/env python

import os

import bioframe as bf
import pandas as pd


def merge_ref_strand_umi(csvfile: os.PathLike) -> pd.DataFrame:
    df = pd.read_csv(
        csvfile, sep="\t", names=["chrom", "start", "end", "umi", "score", "strand"]
    ).drop(columns=["score"])
    df["dsb"] = df.apply(
        lambda row: row["start"] if row["strand"] == "+" else row["end"], axis=1
    )
    df = (
        (
            df.drop(columns=["start", "end"])
            .groupby(["chrom", "umi", "strand"])
            .agg(
                dsb=pd.NamedAgg(column="dsb", aggfunc=lambda x: pd.Series.mode(x)[0]),
                count=pd.NamedAgg(column="dsb", aggfunc="count"),
            )
            .reset_index()
        )
        .rename(columns={"dsb": "start"})
        .assign(end=lambda df: df["start"] + 1)
    )

    columns = df.columns.to_list()
    columns.remove("start")
    columns.insert(1, "start")
    columns.remove("end")
    columns.insert(2, "end")
    df = df.reindex(columns=columns)

    return df


def agg_cluster(df: pd.DataFrame) -> pd.Series:
    count = df["count"].sum()
    vp = ((df["start"] * df["count"]).sum() / count).astype(int)
    return pd.Series(
        {
            "vp": vp,
            "count": count,
        }
    )


def find_cluster(df: pd.DataFrame, cluster_count_thres: int) -> pd.DataFrame:
    df = (
        df.groupby(["chrom", "start", "end"])
        .agg(
            count=pd.NamedAgg(column="count", aggfunc="sum"),
        )
        .reset_index()
        .pipe(
            bf.cluster,
            min_dist=30,
            return_cluster_intervals=False,
        )
    )

    vp = (
        df.groupby(["cluster"])
        .apply(agg_cluster)
        .reset_index()
        .query("count >= @cluster_count_thres")
    ).drop(columns=["count"])

    df = pd.merge(df, vp, on="cluster")
    return df


df = merge_ref_strand_umi("test.csv")
df = find_cluster(df, cluster_count_thres=6)
