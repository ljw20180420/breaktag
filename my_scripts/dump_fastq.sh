#!/bin/bash

function dump_fastq()
{
    # 提取fastq文件
    local sra_acc_list=$1
    local output_dir=$2

    mkdir -p "${cache_dir}/sra"
    while read accsession
    do
        fasterq-dump "${accsession}" --outdir "${output_dir}" --progress
    done < ${sra_acc_list}

    gzip "${cache_dir}/sra/*.fastq"
}
