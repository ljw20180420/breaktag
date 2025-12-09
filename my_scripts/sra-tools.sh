#!/bin/bash

# 切换到当前脚本路径
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# 使用方法：sra-tools.sh sra_acc_list cache_dir

download_sra_until_success() {
    # 第一个参数是SRR的accession编号列表文件
    sra_acc_list=$1

    # 一直下到成功
    until prefetch --max-size u --progress --option-file ${sra_acc_list}
    do
        true
    done
}

validate_sra_until_success() {
    # 第一个参数是SRR的accession编号
    local sra_acc_list=$1
    
    # 一直下到成功
    download_sra_until_success ${sra_acc_list}
    # 一直检查到成功
    until vdb-validate --option-file ${sra_acc_list}
    do
        # 失败了还得重新下载到成功
        download_sra_until_success ${sra_acc_list}
    done
}

sra_acc_list=$1
cache_dir=$2

# 允许缓存
vdb-config -s "cache-enabled=true"

# 下载到user repository
vdb-config  --prefetch-to-user-repo

# 设置user repository路径
vdb-config -s "/repository/user/main/public/root=${cache_dir}"

# 一直检查到成功
validate_sra_until_success ${sra_acc_list} ${cache_dir}
