#!/usr/bin/env python

import os
import subprocess
from io import BytesIO

import pandas as pd


def collect_sgRNA(excelfile: os.PathLike) -> pd.DataFrame:
    def merge_sgRNA(row: pd.Series):
        return "".join(
            row[["Target sequence"] + [f"Unnamed: {i}" for i in range(2, 21)]]
        )

    with pd.ExcelFile(excelfile) as fd:
        print(fd.sheet_names)
        df = (
            pd.read_excel(fd, sheet_name="Table1_sgRNA sequences")
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
    breakpoint()

    return df


df = collect_sgRNA(os.environ["CACHE_DIR"] + "/41587_2024_2238_MOESM3_ESM.xlsx")
