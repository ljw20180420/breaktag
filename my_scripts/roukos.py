#!/usr/bin/env python

import gzip
import os
import pathlib
import re
from collections import OrderedDict

import bioframe as bf
import numpy as np
import pandas as pd
from Bio import Seq

CACHE_DIR = pathlib.Path(os.environ["CACHE_DIR"])


def parse_metadata(matrixfile: os.PathLike) -> pd.DataFrame:

    def get_series_from_title(row: pd.Series) -> str:
        if row["Sample_title"].startswith("U2OS cells"):
            return "48sgrna"
        elif row["Sample_title"].startswith("HepG2 cells") and not row[
            "Sample_title"
        ].startswith("HepG2 cells, LZ3"):
            return "hiplex1"
        elif row["Sample_title"].startswith("HEK293 cells"):
            return "hiplex2"
        elif re.match(r"[FMS]_", row["Sample_title"]):
            return "hiplex3"
        elif row["Sample_title"].startswith("HepG2 cells, LZ3"):
            return "breaktag_1"
        elif re.match(r"(CCR5|DNMT1|EMX|NT|RPL32P3)_", row["Sample_title"]):
            return "breaktag_2"
        elif re.match(r"(PCCF|IDT|PPCF)_", row["Sample_title"]):
            return "cas9_variants"
        else:
            return "other"

    def get_nontarget_from_series(row: pd.Series, df_meta: pd.DataFrame) -> str:
        if row["series"] == "48sgrna":
            nontarget_title = "U2OS cells, Non target bfp gDNA"
        elif row["series"] == "hiplex1":
            nontarget_title = "HepG2 cells, undigested"
        elif row["series"] == "hiplex3":
            nontarget_title = re.sub(
                r"_(REF|ALT)_(pos17|pos18)$", "_NT", row["Sample_title"]
            )
            if nontarget_title not in df_meta["Sample_title"].to_list():
                return "unknown"
        elif row["series"] == "breaktag_1":
            nontarget_title = "HepG2 cells, LZ3, undigested"
        elif row["series"] == "breaktag_2":
            if row["Sample_title"].endswith("_1"):
                nontarget_title = "NT_1"
            else:
                nontarget_title = "NT_2"
        elif row["series"] == "cas9_variants":
            nontarget_title = row["Sample_title"].replace("_pool9Chk", "_NT")
        else:
            return "unknown"

        return df_meta.loc[df_meta["Sample_title"] == nontarget_title, "target"].item()

    with gzip.open(matrixfile, "rb") as fd:
        sample_size = -1
        for line in fd:
            line = line.decode()
            if not line.startswith("!"):
                continue
            sample_size = max(sample_size, len(line.strip().split("\t")) - 1)

    with gzip.open(matrixfile, "rb") as fd:
        dc = {}
        for line in fd:
            line = line.decode()
            if not line.startswith("!"):
                continue
            fields = line.strip().split("\t")
            head, fields = fields[0], [field.strip('"') for field in fields[1:]]
            if len(fields) < sample_size:
                continue
            dc[head.lstrip("!")] = fields

    return (
        pd.DataFrame(dc)
        .query('Sample_supplementary_file_1.str.endswith(".bed.gz")')
        .assign(
            series=lambda df: df.apply(get_series_from_title, axis=1),
            target=lambda df: df["Sample_supplementary_file_1"].map(os.path.basename),
            nontarget=lambda df: df.apply(
                lambda row, df_meta=df: get_nontarget_from_series(row, df_meta), axis=1
            ),
        )
        .query("target != nontarget")
        .reset_index(drop=True)
    )


def collect_sgRNA(excelfile: os.PathLike) -> pd.DataFrame:
    def merge_sgRNA(row: pd.Series):
        return "".join(
            row[["Target sequence"] + [f"Unnamed: {i}" for i in range(2, 21)]]
        )

    with pd.ExcelFile(excelfile) as fd:
        df_sgRNA = (
            pd.read_excel(fd, sheet_name="Table1_sgRNA sequences")
            .drop_duplicates(subset=["sgRNA"], keep="first", ignore_index=True)
            .assign(seq=lambda df: df.apply(merge_sgRNA, axis=1))
            .drop(columns=["Target sequence"] + [f"Unnamed: {i}" for i in range(2, 21)])
        )

    return df_sgRNA


def get_sgRNA(row: pd.Series, df_sgRNA: pd.DataFrame):
    if row["series"] == "48sgrna":
        return df_sgRNA.loc[
            "U2OS cells, " + df_sgRNA["sgRNA"] + " gDNA" == row["Sample_title"],
            "seq",
        ].item()

    if row["series"] == "hiplex1":
        Dataset_value = "Pool{}_Hiplex1".format(
            re.search(r"HepG2 cells, Pool of sgRNAs no\. (\w+)$", row["Sample_title"])
            .group(1)
            .capitalize()
        )
    elif row["series"] == "hiplex2":
        Dataset_value = "Pool{}_Hiplex2".format(
            re.search(r"HEK293 cells, Pool of sgRNAs no\. (\w+)$", row["Sample_title"])
            .group(1)
            .capitalize()
        )
    elif row["series"] == "hiplex3":
        Dataset_value = "REF_{}_Hiplex3".format(
            re.search(r"_(pos17|pos18)$", row["Sample_title"]).group(1)
        )
    else:
        return "unknown"

    return ":".join(df_sgRNA.loc[df_sgRNA["Dataset"] == Dataset_value, "seq"])


def load_fasta(genomefile: os.PathLike) -> OrderedDict:
    genomefile = pathlib.Path(os.fspath(genomefile))
    refseq2ucsc = pd.read_json(
        genomefile.parent / "sequence_report.jsonl", lines=True
    ).set_index("refseqAccession")["ucscStyleName"]
    genome = bf.load_fasta(os.fspath(genomefile))
    genome = OrderedDict(
        [
            (refseq2ucsc[key], val)
            for key, val in genome.items()
            if re.search(r"^chr(\d{1,2}|X|Y|M)$", refseq2ucsc[key])
        ]
    )
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

    sense_mismatches = (sgRNA_lib == np.array(list(sense_sgRNA))[None, :]).sum(axis=1)
    sense_min_idx = sense_mismatches.argmin()
    antisense_mismatches = (sgRNA_lib == np.array(list(antisense_sgRNA))[None, :]).sum(
        axis=1
    )
    antisense_min_idx = antisense_mismatches.argmin()

    sense = (
        "+"
        if sense_mismatches[sense_min_idx] <= antisense_mismatches[antisense_min_idx]
        else "-"
    )

    return pd.Series(
        {
            "sense": sense,
            "sgRNA": sense_sgRNA if sense == "+" else antisense_sgRNA,
            "PAM": (
                chromsome.ff.fetch(
                    chromsome.ref, row["end"] + 3, row["end"] + 6
                ).upper()
                if sense == "+"
                else str(
                    Seq.Seq(
                        chromsome.ff.fetch(
                            chromsome.ref, row["end"] - 6, row["end"] - 3
                        )
                    )
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
        }
    )


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
            pd.read_csv(
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
        pd.concat(
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

    total_breaks = pd.concat(
        [
            breaks["target"].assign(
                start=lambda df: df["end"], end=lambda df: df["end"] + 1, tnt="target"
            ),
            breaks["nontarget"].assign(
                start=lambda df: df["end"],
                end=lambda df: df["end"] + 1,
                tnt="nontarget",
            ),
        ]
    )

    offtarget = (
        bf.overlap(offtarget, total_breaks)
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
        pd.pivot_table(
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


df_meta = (
    pd.concat(
        [
            parse_metadata(CACHE_DIR / "GSE223772-GPL18573_series_matrix.txt.gz"),
            parse_metadata(CACHE_DIR / "GSE223772-GPL24676_series_matrix.txt.gz"),
        ]
    )
    .reset_index(drop=True)
    .replace("", np.nan)
    .dropna(axis=1, how="any")
)
df_meta.to_csv(CACHE_DIR / "meta.csv", index=False)

df_sgRNA = collect_sgRNA(CACHE_DIR / "41587_2024_2238_MOESM3_ESM.xlsx")
df_sgRNA.to_csv(CACHE_DIR / "sgRNA.csv", index=False)

df_meta = df_meta.assign(
    sgRNA_lib=lambda df: df.apply(
        lambda row, df_sgRNA=df_sgRNA: get_sgRNA(row, df_sgRNA), axis=1
    )
)

genome = load_fasta(
    CACHE_DIR
    / "human_dataset"
    / "ncbi_dataset"
    / "data"
    / "GCF_000001405.40"
    / "GCF_000001405.40_GRCh38.p14_genomic.fna"
)


# # test code
# target = "/home/ljw/new_fold/projects/roukoslab/breakinspectoR/.conda/lib/R/library/breakinspectoR/extdata/vegfa.chr6.bed.gz"
# nontarget = "/home/ljw/new_fold/projects/roukoslab/breakinspectoR/.conda/lib/R/library/breakinspectoR/extdata/nontarget.chr6.bed.gz"
# sgRNA_lib = np.array([list("GACCCCCTCCACCCCGCCTC")])
# min_breaks = 2
# ext = 17
# max_mismatch = 7
# offtarget = analyze_bed(
#     target, nontarget, sgRNA_lib, min_breaks, max_mismatch, genome, ext
# )

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
    sgRNA_lib = np.array([list(sgRNA) for sgRNA in sgRNA_lib.split(":")])
    offtarget = analyze_bed(
        target, nontarget, sgRNA_lib, min_breaks, max_mismatch, genome, ext
    )
    offtarget.to_csv(CACHE_DIR / "result" / f"{sample_title}.{series}.csv", index=False)
