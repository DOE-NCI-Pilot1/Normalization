#!/bin/bash

# bash scripts/run.bash
datapath=data/Dec2019/combined_rnaseq_data_lincs1000
# src="ccle nci60 ncipdm nova"
src="ccle nci60 ncipdm"
outpath=out/Dec2019/combined_rnaseq_data_lincs1000

python src/normalize_rnaseq.py \
    --datapath $datapath \
    --outpath $outpath \
    --src $src
