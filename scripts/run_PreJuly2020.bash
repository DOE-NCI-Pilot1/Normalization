#!/bin/bash

# Example:
# bash scripts/run.bash

dname=PreJuly2020/combined_rnaseq_data_lincs1000
# dname=PreJuly2020/combined_rnaseq_data

datapath=data/$dname
outpath=out/$dname

src="ccle nci60 ncipdm gdc"
# src="ccle nci60 ncipdm"

python src/normalize_rnaseq.py \
    --datapath $datapath \
    --outpath $outpath \
    --src $src
