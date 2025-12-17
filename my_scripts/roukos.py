#!/usr/bin/env python

import gzip
import os
import re
import subprocess
from io import BytesIO

import bioframe as bf
import numpy as np
import pandas as pd


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
    if row["series"] == "hiplex1":
        nontarget_title = "HepG2 cells, undigested"
    elif row["series"] == "hiplex3":
        nontarget_title = re.sub(
            r"_(REF|ALT)_(pos17|pos18)$", "_NT", row["Sample_title"]
        )
        if nontarget_title not in df_meta["Sample_title"]:
            nontarget_title = "unknown"
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
        nontarget_title = "unknown"

    if nontarget_title == "unknown":
        return "unknown"

    return df_meta.loc[df_meta["Sample_title"] == nontarget_title, "target"].item()


def parse_metadata(matrixfile: os.PathLike) -> pd.DataFrame:
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


def analyze_bed(target: os.PathLike, nontarget: os.PathLike) -> pd.DataFrame:
    df_target = pd.read_csv(
        target, sep="\t", names=["chrom", "start", "end", "name", "score", "strand"]
    )
    df_nontarget = pd.read_csv(
        nontarget, sep="\t", names=["chrom", "start", "end", "name", "score", "strand"]
    )

    breakpoint()


def collect_sgRNA(excelfile: os.PathLike) -> pd.DataFrame:
    def merge_sgRNA(row: pd.Series):
        return "".join(
            row[["Target sequence"] + [f"Unnamed: {i}" for i in range(2, 21)]]
        )

    with pd.ExcelFile(excelfile) as fd:
        print(fd.sheet_names)
        df = (
            pd.read_excel(fd, sheet_name="Table1_sgRNA sequences")
            .drop_duplicates(subset=["sgRNA"], keep="first", ignore_index=True)
            .assign(seq=lambda df: df.apply(merge_sgRNA, axis=1))
            .drop(columns=["Target sequence"] + [f"Unnamed: {i}" for i in range(2, 21)])
        )

    # 基因组匹配
    sgRNA_fasta_file = os.environ["CACHE_DIR"] + "/result/sgRNA.fasta"
    with open(sgRNA_fasta_file, "w") as fd:
        for name, seq in zip(df["sgRNA"], df["seq"]):
            fd.write(f">{name}\n{seq}\n")
    genome = (
        os.environ["CACHE_DIR"]
        + "/human_dataset/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"
    )
    result = subprocess.run(
        args=f"bwa mem -a -v 0 -T 20 {genome} {sgRNA_fasta_file} | bedtools bamtobed -i -",
        shell=True,
        capture_output=True,
    )

    # 找到vp
    seq_report = pd.read_json(
        os.environ["CACHE_DIR"]
        + "/human_dataset/ncbi_dataset/data/GCF_000001405.40/sequence_report.jsonl",
        lines=True,
    )
    vp_df = (
        pd.read_csv(
            BytesIO(result.stdout),
            sep="\t",
            names=["chrom", "start", "end", "sgRNA", "mapQ", "strand"],
        )
        .assign(
            chrom=lambda df: df["chrom"].map(
                seq_report.set_index("refseqAccession")["ucscStyleName"]
            )
        )
        .query('chrom.str.match(r"^chr([0-9]+|X|Y|M)$")')
        .reset_index(drop=True)
    )

    def pam_range(row: pd.Series) -> pd.Series:
        return pd.Series(
            {
                "chrom": row["chrom"],
                "start": row["end"] if row["strand"] == "+" else row["start"] - 6,
                "end": row["start"] if row["strand"] == "-" else row["end"] + 6,
                "name": row["sgRNA"],
                "score": row["mapQ"],
                "strand": row["strand"],
            }
        )

    pam_bed_file = os.environ["CACHE_DIR"] + "/result/pam.bed"
    vp_df.apply(pam_range, axis=1).assign(
        chrom=lambda df: df["chrom"].map(
            seq_report.set_index("ucscStyleName")["refseqAccession"]
        )
    ).to_csv(pam_bed_file, index=False, header=False, sep="\t")

    result = subprocess.run(
        args=f"bedtools getfasta -fi {genome} -bed {pam_bed_file} -bedOut -s",
        shell=True,
        capture_output=True,
    )

    vp_df["pam"] = pd.read_csv(
        BytesIO(result.stdout),
        sep="\t",
        names=["chrom", "start", "end", "sgRNA", "mapQ", "strand", "pam"],
    )["pam"].str.upper()

    vp_df.query('pam.str.match(r".GG")')

    return df


df_meta = (
    pd.concat(
        [
            parse_metadata(
                os.environ["CACHE_DIR"] + "/GSE223772-GPL18573_series_matrix.txt.gz"
            ),
            parse_metadata(
                os.environ["CACHE_DIR"] + "/GSE223772-GPL24676_series_matrix.txt.gz"
            ),
        ]
    )
    .reset_index(drop=True)
    .replace("", np.nan)
    .dropna(axis=1, how="any")
)
df_meta.to_csv(os.environ["CACHE_DIR"] + "/meta.csv", index=False)

for sample_title, series, target, nontarget in zip(
    df_meta["Sample_title"], df_meta["series"], df_meta["target"], df_meta["nontarget"]
):
    if not target.endswith(".bed.gz") or nontarget == "unknown":
        continue

    analyze_bed(target, nontarget)

# 48sgrna, hiplex1, hiplex2, hiplex3, breaktag_1, breaktag_2, cas9_variants
# collect_bedgz("48sgrna", df_meta)
# df = collect_sgRNA(os.environ["CACHE_DIR"] + "/41587_2024_2238_MOESM3_ESM.xlsx")
