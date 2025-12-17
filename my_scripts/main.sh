#!/bin/bash

# 切换到当前脚本路径
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) 

. ./dump_fastq.sh
. ./download_genome.sh

sra_acc_list="${1:-"SRA_ACCESSION_LIST_FILE.txt"}"
gse_acc_list="${2:-"my_scripts/GSE.txt"}"
cache_dir="${3:-"${CACHE_DIR}"}"

# # 下载sra提取fastq并压缩.
# dump_fastq ${sra_acc_list} ${cache_dir}/fastq

# # 下载基因组并解压.
# download_accession "GCF_000001405.40" "${cache_dir}/human_dataset.zip"

# # bwa index基因组.
# bwa index "${cache_dir}/human_dataset/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"

# # 下载suppl文件
# wget -c -O ${cache_dir}/breaktag_raw_data/GSE223772_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE223772&format=file"

# # 下载metadata
# wget -c -P ${cache_dir}/ https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223772/matrix/GSE223772-GPL18573_series_matrix.txt.gz
# wget -c -P ${cache_dir}/ https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223772/matrix/GSE223772-GPL24676_series_matrix.txt.gz

# # 运行pipeline.
# mkdir -p "${cache_dir}/test"
# cp ../pipelines/breaktag/targets.txt "${cache_dir}/test"
# bpipe run ../pipelines/breaktag/breaktag.pipeline.groovy "${cache_dir}/fastq/SRR23241572.fastq.gz"
