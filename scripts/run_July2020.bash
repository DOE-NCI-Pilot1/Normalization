#!/bin/bash

# Example:
# bash scripts/run.bash

dname=July2020/combined_rnaseq_data_lincs1000
# dname=July2020/combined_rnaseq_data

datapath=data/$dname
outpath=out/$dname

src="ccle nci60 ncipdm nova gdc"
# src="ccle nci60 ncipdm nova"
# src="ccle nci60 ncipdm gdc"
# src="ccle nci60 ncipdm"

python src/normalize_rnaseq.py \
    --datapath $datapath \
    --outpath $outpath \
    --src $src
