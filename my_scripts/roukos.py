#!/usr/bin/env python

import os
import pathlib
import re
from collections import OrderedDict

import bioframe as bf
import numpy as np
import pandas as pd
from Bio import Seq

CACHE_DIR = pathlib.Path(os.environ["CACHE_DIR"])


def load_fasta(genomefile: os.PathLike) -> OrderedDict:
    genomefile = pathlib.Path(os.fspath(genomefile))
    refseq2ucsc = pd.read_json(
        genomefile.parent / "sequence_report.jsonl", lines=True
    ).set_index("refseqAccession")["ucscStyleName"]
    genome = bf.load_fasta(os.fspath(genomefile))
    genome = OrderedDict([
        (refseq2ucsc[key], val)
        for key, val in genome.items()
        if re.search(r"^chr(\d{1,2}|X|Y|M)$", refseq2ucsc[key])
    ])
    return genome


def filter_boundary_entry(row: pd.Series, genome: OrderedDict):
    chromsome = genome[row["chrom"]]
    chromsome_length = chromsome.ff.get_reference_length(chromsome.ref)
    return row["end"] - 17 >= 0 and row["end"] + 17 < chromsome_length


def get_sense_sgRNA_pam(
    row: pd.Series, genome: OrderedDict, sgRNA_lib: np.ndarray
) -> pd.Series:
    chromsome = genome[row["chrom"]]
    sense_sgRNA = chromsome.ff.fetch(
        chromsome.ref, row["end"] - 17, row["end"] + 3
    ).upper()
    antisense_sgRNA = str(
        Seq.Seq(
            chromsome.ff.fetch(chromsome.ref, row["end"] - 3, row["end"] + 17)
        ).reverse_complement()
    ).upper()

    sense_mismatches = (sgRNA_lib != np.array(list(sense_sgRNA))[None, :]).sum(axis=1)
    sense_min_idx = sense_mismatches.argmin()
    antisense_mismatches = (sgRNA_lib != np.array(list(antisense_sgRNA))[None, :]).sum(
        axis=1
    )
    antisense_min_idx = antisense_mismatches.argmin()

    sense = (
        "+"
        if sense_mismatches[sense_min_idx] <= antisense_mismatches[antisense_min_idx]
        else "-"
    )

    return pd.Series({
        "sense": sense,
        "sgRNA": sense_sgRNA if sense == "+" else antisense_sgRNA,
        "PAM": (
            chromsome.ff.fetch(chromsome.ref, row["end"] + 3, row["end"] + 6).upper()
            if sense == "+"
            else str(
                Seq
                .Seq(chromsome.ff.fetch(chromsome.ref, row["end"] - 6, row["end"] - 3))
                .reverse_complement()
                .upper()
            )
        ),
        "sgRNA_lib": "".join(
            sgRNA_lib[sense_min_idx if sense == "+" else antisense_min_idx]
        ),
        "sgRNA_mismatch": (
            sense_mismatches[sense_min_idx]
            if sense == "+"
            else antisense_mismatches[antisense_min_idx]
        ),
    })


def analyze_bed(
    target: os.PathLike,
    nontarget: os.PathLike,
    sgRNA_lib: np.ndarray,
    min_breaks: int,
    max_mismatch: int,
    genome: OrderedDict,
    ext: int,
) -> pd.DataFrame:
    breaks = {}
    for name, bedfile in zip(["target", "nontarget"], [target, nontarget]):
        breaks[name] = (
            pd
            .read_csv(
                bedfile,
                sep="\t",
                names=["chrom", "start", "end", "name", "score", "strand"],
            )
            .query('chrom.str.match(r"^chr([0-9]{1,2}|X|Y|M)$")')
            .groupby(["chrom", "start", "end", "strand"])
            .agg(score=pd.NamedAgg(column="score", aggfunc="sum"))
            .reset_index()
        )

    offtarget = (
        breaks["target"]
        .groupby(["chrom", "start", "end"])
        .agg(score=pd.NamedAgg(column="score", aggfunc="sum"))
        .query("score >= @min_breaks")
        .reset_index()
        .drop(columns=["score"])
    )
    offtarget = offtarget.loc[
        offtarget.apply(
            lambda row, genome=genome: filter_boundary_entry(row, genome), axis=1
        )
    ]
    offtarget = (
        pd
        .concat(
            [
                offtarget,
                offtarget.apply(
                    lambda row, genome=genome, sgRNA_lib=sgRNA_lib: get_sense_sgRNA_pam(
                        row, genome, sgRNA_lib
                    ),
                    axis=1,
                ),
            ],
            axis=1,
        )
        .query("sgRNA_mismatch <= @max_mismatch")
        .assign(
            cut=lambda df: df["end"],
            start=lambda df: df["cut"] - ext,
            end=lambda df: df["cut"] + ext + 1,
        )
    )

    total_breaks = pd.concat([
        breaks["target"].assign(
            start=lambda df: df["end"],
            end=lambda df: df["end"] + 1,
            tnt="target",
        ),
        breaks["nontarget"].assign(
            start=lambda df: df["end"],
            end=lambda df: df["end"] + 1,
            tnt="nontarget",
        ),
    ])

    offtarget = (
        bf
        .overlap(offtarget, total_breaks)
        .assign(
            break_pos=lambda df: df.apply(
                lambda row: (
                    row["start_"] - row["cut"]
                    if row["sense"] == "+"
                    else row["cut"] - row["start_"]
                ),
                axis=1,
            ).astype(
                pd.CategoricalDtype(categories=range(-ext, ext + 1), ordered=True)
            ),
        )
        .rename(columns={"strand_": "strand", "score_": "score", "tnt_": "tnt"})
        .drop(columns=["chrom_", "start_", "end_"])
    )

    index_columns = [
        column
        for column in offtarget.columns.to_list()
        if column not in ["score", "break_pos"]
    ]

    offtarget = (
        pd
        .pivot_table(
            offtarget,
            values="score",
            index=index_columns,
            columns=["break_pos"],
            aggfunc="sum",
            fill_value=0,
            observed=True,
        )
        .reset_index()
        .reindex(columns=index_columns + list(range(-ext, ext + 1)))
        .fillna(0)
        .astype({pos: int for pos in range(-ext, ext + 1)})
    )

    return offtarget


if __name__ == "__main__":
    df_meta = pd.read_csv(CACHE_DIR / "meta.csv", header=0)

    genome = load_fasta(
        CACHE_DIR
        / "human_dataset"
        / "ncbi_dataset"
        / "data"
        / "GCF_000001405.40"
        / "GCF_000001405.40_GRCh38.p14_genomic.fna"
    )

    min_breaks = 2
    ext = 17
    max_mismatch = 7
    for sample_title, series, target, nontarget, sgRNA_lib in zip(
        df_meta["Sample_title"],
        df_meta["series"],
        df_meta["target"],
        df_meta["nontarget"],
        df_meta["sgRNA_lib"],
    ):
        if nontarget == "unknown" or sgRNA_lib == "unknown":
            continue

        target = CACHE_DIR / "breaktag_raw_data" / target
        nontarget = CACHE_DIR / "breaktag_raw_data" / nontarget
        breakpoint()
        sgRNA_lib = np.array([list(sgRNA) for sgRNA in sgRNA_lib.split(":")])
        offtarget = analyze_bed(
            target, nontarget, sgRNA_lib, min_breaks, max_mismatch, genome, ext
        )
        offtarget.to_csv(
            CACHE_DIR / "result" / f"{sample_title}.{series}.csv", index=False
        )
