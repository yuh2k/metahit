#!/usr/bin/env bash

# 使用示例：
# ./bin_refinement.sh --fasta reference.fasta --contig_file contig_info.p --contact_matrix_file contact_matrix.npz --output output_dir

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 [options]"
    exit 1
fi

# 解析命令行参数
./bin/metahit-scripts/bin_refinement.py "$@"