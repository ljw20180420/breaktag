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
# while read accsession
# do
#     bedfile=$(
#         curl https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6995nnn/${accsession}/suppl/ |
#         sed '$a </body></html>' |
#         xq -m |
#         xq -x /html/body/pre/a |
#         grep ${accsession}
#     )
    
#     wget -P ${cache_dir}/breaktag_raw_data/ https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6995nnn/${accsession}/suppl/${bedfile}
# done < GSE.txt

# # 运行pipeline.
# mkdir -p "${cache_dir}/test"
# cp ../pipelines/breaktag/targets.txt "${cache_dir}/test"
# bpipe run ../pipelines/breaktag/breaktag.pipeline.groovy "${cache_dir}/fastq/SRR23241572.fastq.gz"
