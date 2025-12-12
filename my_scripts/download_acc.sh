#!/bin/bash

prefetch_acc_until_success() {
    # accession编号列表文件
    local acc_list=$1
    local output_dir=$2

    # 一直下到成功
    until prefetch --max-size u --progress --output-directory ${output_dir} --option-file ${acc_list}
    do
        true
    done
}

validate_acc_until_success() {
    # accession编号列表文件
    local acc_list=$1
    local output_dir=$2

    # 一直下到成功
    download_acc_until_success ${acc_list} ${output_dir}
    # 一直检查到成功
    until vdb-validate --option-file ${acc_list} ${output_dir}
    do
        # 失败了还得重新下载到成功
        download_acc_until_success ${acc_list} ${output_dir}
    done
}

download_acc() {
    local acc_list=$1
    local output_dir=$2

    # 一直检查到成功
    validate_acc_until_success ${acc_list} ${output_dir}
}
