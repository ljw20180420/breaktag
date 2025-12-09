#!/bin/bash

# 切换到当前脚本路径
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# 使用方法：fasterq_dump.sh sra_acc_list cache_dir

# 提取fastq文件
sra_acc_list=$1
cache_dir=$2
while read accsession
do
    fasterq-dump $accsession --outdir ${cache_dir}/sra --progress
done < ${sra_acc_list}
