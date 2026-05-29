#!/usr/bin/env python

import gzip
import os
import pathlib
import re

import numpy as np
import pandas as pd

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
        pd
        .DataFrame(dc)
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


def get_sgRNA(row: pd.Series, df_sgRNA: pd.DataFrame):
    if row["series"] == "48sgrna":
        return df_sgRNA.loc[
            "U2OS cells, " + df_sgRNA["sgRNA"] + " gDNA" == row["Sample_title"],
            "seq",
        ].item()

    if row["series"] == "hiplex1":
        Dataset_value = "Pool{}_Hiplex1".format(
            re
            .search(r"HepG2 cells, Pool of sgRNAs no\. (\w+)$", row["Sample_title"])
            .group(1)
            .capitalize()
        )
    elif row["series"] == "hiplex2":
        Dataset_value = "Pool{}_Hiplex2".format(
            re
            .search(r"HEK293 cells, Pool of sgRNAs no\. (\w+)$", row["Sample_title"])
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


def collect_sgRNA(excelfile: os.PathLike) -> pd.DataFrame:
    def merge_sgRNA(row: pd.Series):
        return "".join(
            row[["Target sequence"] + [f"Unnamed: {i}" for i in range(2, 21)]]
        )

    with pd.ExcelFile(excelfile) as fd:
        df_sgRNA = (
            pd
            .read_excel(fd, sheet_name="Table1_sgRNA sequences")
            .drop_duplicates(subset=["sgRNA"], keep="first", ignore_index=True)
            .assign(seq=lambda df: df.apply(merge_sgRNA, axis=1))
            .drop(columns=["Target sequence"] + [f"Unnamed: {i}" for i in range(2, 21)])
        )

    return df_sgRNA


if __name__ == "__main__":
    df_meta = (
        pd
        .concat([
            parse_metadata(CACHE_DIR / "GSE223772-GPL18573_series_matrix.txt.gz"),
            parse_metadata(CACHE_DIR / "GSE223772-GPL24676_series_matrix.txt.gz"),
        ])
        .reset_index(drop=True)
        .replace("", np.nan)
        .dropna(axis=1, how="any")
    )

    df_sgRNA = collect_sgRNA(CACHE_DIR / "41587_2024_2238_MOESM3_ESM.xlsx")
    df_sgRNA.to_csv(CACHE_DIR / "sgRNA.csv", index=False)

    df_meta = df_meta.assign(
        sgRNA_lib=lambda df: df.apply(
            lambda row, df_sgRNA=df_sgRNA: get_sgRNA(row, df_sgRNA), axis=1
        )
    )
    df_meta.to_csv(CACHE_DIR / "meta.csv", index=False)
