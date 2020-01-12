import warnings
warnings.filterwarnings('ignore')

import os
import sys
import argparse
from pathlib import Path
from time import time
from glob import glob
from pprint import pprint

import sklearn
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt

# Default settings
default_datapath = Path('../data/combined_rnaseq_data_lincs1000')
default_sources = ['ccle', 'nci60', 'ncipdm', 'gdc']
src_name_map = {'ccle': 'CCLE', 'nci60': 'NCI-60', 'ncipdm': 'NCI-PDM', 'gdc': 'GDC'}
dpi = 100

# Utils
from utils import scale_rna, plot_pca, load_rna, combat_

# File path
filepath = Path(__file__).resolve().parent



def src_from_cell_col(cell_name_col, verbose=False):
    """ Takes column that contains sample names, extract the source name,
    and returns the arr of source names. Prints value_counts if verbose=True.
    Args:
        cell_name_col : the actual array of sample names
    """
    src_names_arr = cell_name_col.map(lambda x: x.split('.')[0].lower())
    src_names_arr.name = 'source'
    if verbose:
        print(src_names_arr.value_counts())
    return src_names_arr



def parse_args(args):
    parser = argparse.ArgumentParser(description="Normalize RNA-Seq (batch-effect removal).")
    parser.add_argument('-dp', '--datapath', default=default_datapath, type=str, help='Full path to RNA-Seq data (default: None).')
    parser.add_argument('--src', nargs='+', default=default_sources, type=str, help='List of sources to remove batch effect (default: None).')
    parser.add_argument('--fea_prfx', default=None, type=str, help='Prefix to add to each feature (default: None).')
    parser.add_argument('-o', '--outpath', default='./out_figs', type=str, help="Output path to store figures (default: './out_figs').")

    args = parser.parse_args(args)
    return args


def run(args):
    fea_start_id = 1 # col index of the first gene feature
    outpath = Path(args['outpath'])
    fea_prfx = args['fea_prfx']
    os.makedirs(outpath, exist_ok=True)

    # ---------------------------------------------------------
    #   Load data
    # ---------------------------------------------------------
    print('\nLoad RNA-Seq ... {}'.format( args['datapath']) )
    rna = pd.read_csv(Path(args['datapath']), sep='\t', na_values=['na', '-', ''], warn_bad_lines=True)
    if fea_prfx is not None:
        rna = rna.rename(columns={c: fea_prfx+c for c in rna.columns[1:] if fea_prfx not in c}) # prefix rna gene names
    rna.rename(columns={'Sample': 'CELL'}, inplace=True)

    print(f'rna.shape {rna.shape}')
    # src_from_cell_col(rna['CELL'], verbose=True);


    # ---------------------------------------------------------
    #   Extract sources which contain original RNA-Seq
    # ---------------------------------------------------------
    src_names_arr = src_from_cell_col(rna['CELL'], verbose=True)
    rna = rna.loc[ src_names_arr.isin(args['src']) ]
    print(f'rna.shape {rna.shape}')
    # src_from_cell_col(rna['CELL'], verbose=True);


    # ---------------------------------------------------------
    #   Raw data
    # ---------------------------------------------------------
    print('\nRaw Data ...')
    src_names_arr = src_from_cell_col(rna['CELL'], verbose=True)

    plot_pca(rna.iloc[:, fea_start_id:], components = [1, 2], figsize=(8,7),
             # color_vector = rna['CELL'].map(lambda x: src_name_map[x.split('.')[0].lower()]),
             color_vector = src_names_arr,
             scale=True, title='PCA (raw data)');
    plt.savefig(outpath/'pca_raw.png', dpi=dpi)


    # ---------------------------------------------------------
    #   Source scaling
    # ---------------------------------------------------------
    print('\nSource scaling ...')

    # Source scale
    rna_src = scale_rna(rna, fea_start_id=fea_start_id, per_source=True)
    print(rna_src.shape)

    plot_pca(rna_src.iloc[:, fea_start_id:], components = [1, 2], figsize=(8,7),
             # color_vector = rna['CELL'].map(lambda x: src_name_map[x.split('.')[0].lower()]),
             color_vector = src_names_arr,
             scale=True, title='PCA (source-scaled data)');
    plt.savefig(outpath/'pca_src_scale.png', dpi=dpi)


    # ---------------------------------------------------------
    #   Combat scaling
    # ---------------------------------------------------------
    print('\nComBat scaling ...')

    # Data for combat
    meta = src_from_cell_col(rna['CELL'], verbose=False).to_frame()
    meta['CELL'] = rna['CELL']

    # Combat scale
    rna_combat = combat_(rna, meta, sample_col_name='CELL', batch_col_name='source')

    plot_pca(rna_combat.iloc[:, fea_start_id:], components = [1, 2], figsize=(8,7),
             # color_vector = rna['CELL'].map(lambda x: src_name_map[x.split('.')[0].lower()]),
             color_vector = src_names_arr,
             scale=True, title='PCA (combat-scaled data)');
    plt.savefig(outpath/'pca_combat.png', dpi=dpi)


def main(args):
    args = parse_args(args)
    args = vars(args)
    run(args)
    print('Done.')


if __name__ == '__main__':
    main(sys.argv[1:])


