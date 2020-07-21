import warnings
warnings.filterwarnings('ignore')

import os
import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt

from sklearn.impute import SimpleImputer

# Default settings
default_datapath = Path('./data/Dec2019/combined_rnaseq_data_lincs1000')
# default_sources = ['ccle', 'nci60', 'ncipdm', 'gdc', 'nova']
default_sources = ['ccle', 'nci60', 'ncipdm', 'nova']
src_name_map = {'ccle': 'CCLE', 'nci60': 'NCI-60', 'ncipdm': 'NCI-PDM',
                'gdc': 'GDC', 'nova': 'Novartis'}
dpi = 100

# Utils
from utils import scale_rna, plot_pca, load_rna, combat_

# File path
filepath = Path(__file__).resolve().parent

OUTPATH = Path(filepath, 'out')


def parse_args(args):
    parser = argparse.ArgumentParser(description='Normalize RNA-Seq (batch-effect removal).')

    parser.add_argument('-dp', '--datapath',
                        default=default_datapath,
                        type=str,
                        help=f'Full path to RNA-Seq data (default: {default_datapath}).')
    parser.add_argument('-o', '--outpath',
                        default=OUTPATH,
                        type=str,
                        help=f'Output path to store figures (default: {OUTPATH}).')
    parser.add_argument('--src',
                        nargs='+',
                        default=default_sources,
                        type=str,
                        help=f'List of sources to remove batch effect (default: {default_sources}).')
    parser.add_argument('--fea_prfx',
                        default=None,
                        type=str,
                        help='Prefix to add to each feature (default: None).')
    parser.add_argument('--include_gdc',
                        action='store_true',
                        help='Whether to include GDC (TCGA) data (default: False).')

    args = parser.parse_args(args)
    return args


def src_from_cell_col(cell_name_col, verbose=False):
    """ Extract the source name from sample names. 
    Args:
        cell_name_col: the actual array of sample names
        verbose: print value_counts if True
    """
    src_names_arr = cell_name_col.map(lambda x: x.split('.')[0].lower())
    src_names_arr.name = 'source'
    if verbose:
        print(src_names_arr.value_counts())
    return src_names_arr


def map_rna_profiles(rna, cl_map, sample_col_name='Sample'):
    """
    Use a mapping file (cl_map) to copy ccle and nci60 gene expression to gdsc, gcsi, and ctrp.
    Args:
        rna: df that contains only data sources that have their original rna expression data
        cl_map: df with columns 'from_cell' and 'to_cell' 
    Returns:
        df: contains the merged df which contains the original and the mapped samples
    """
    # Merge in order to copy rnaseq profiles
    cells = cl_map.merge(rna, left_on='from_cell', right_on=sample_col_name, how='inner')
    
    # Drop and rename columns
    cells = cells.drop(columns=['from_cell', sample_col_name])
    cells = cells.rename(columns={'to_cell': sample_col_name})

    # Concat 'rna' (rna that contains original rna expression) and 'cells' replicated expression
    df = pd.concat([rna, cells], axis=0).sort_values(sample_col_name).reset_index(drop=True)
    return df


def run(args):
    import ipdb; ipdb.set_trace(context=5)
    fea_start_id = 1  # col index of the first gene feature
    datapath = Path(args.datapath)
    outpath = Path(args.outpath)
    fea_prfx = args.fea_prfx
    sample_col_name = 'Sample'
    # sample_col_name = 'CELL'
    os.makedirs(outpath, exist_ok=True)

    # ---------------------------------------------------------
    #   Load data
    # ---------------------------------------------------------
    print('\nLoad RNA-Seq ... {}'.format(datapath))
    rna = pd.read_csv(datapath, sep='\t', na_values=['na', '-', ''])
    if fea_prfx is not None:
        rna = rna.rename(columns={c: fea_prfx+c for c in rna.columns[1:] if fea_prfx not in c})  # prfx gene names
    # rna.rename(columns={'Sample': 'CELL'}, inplace=True)
    print(f'rna.shape {rna.shape}')

    # Load cell line mapping
    cl_mapping = pd.read_csv('data/cl_mapping', sep='\t', names=['from_cell', 'to_cell'], header=0)

    # Keep or drop GDC (TCGA)
    if args.include_gdc is False:
        aa = src_from_cell_col(rna[sample_col_name], verbose=True);
        idx = (aa != 'gdc')  # non-gdc indices
        rna = rna[idx].reset_index(drop=True)
        src_from_cell_col(rna[sample_col_name], verbose=True);

    # ---------------------------------------------------------
    #   Extract sources which contain original RNA-Seq
    # ---------------------------------------------------------
    rna_copy = rna.copy()
    src_names_arr = src_from_cell_col(rna[sample_col_name], verbose=True)
    rna = rna_copy.loc[ src_names_arr.isin(args.src) ].reset_index(drop=True)
    print(f'rna.shape {rna.shape}')
    src_from_cell_col(rna[sample_col_name], verbose=True);

    # ---------------------------------------------------------
    #   Impute missing values
    # ---------------------------------------------------------
    print('\nTotal columns with NaNs: {}'.format( sum(rna.iloc[:, 1:].isna().sum() > 0) ))

    # Dump the count of NANs in each column
    aa = rna.isna().sum(axis=0).reset_index()
    aa = aa.rename(columns={'index': 'col', 0: 'count'})
    aa = aa.sort_values('count', ascending=False).reset_index(drop=True)
    aa.to_csv(outpath/'nan_count_per_col.csv', index=False)

    # Data imputation
    # TODO: should we consider a different approach than the mean??
    if sum( rna.isna().sum() > 0 ):
        imputer = SimpleImputer(missing_values=np.nan, strategy='mean')
        rna.iloc[:, 1:] = imputer.fit_transform(rna.iloc[:, 1:].values)
        print('Imputed missing values.')
        print('Total columns with NaNs: {}'.format( sum(rna.iloc[:, 1:].isna().sum() > 0) ))

    # ---------------------------------------------------------
    #   Raw data
    # ---------------------------------------------------------
    print('\nRaw Data ...')
    src_names_arr = src_from_cell_col(rna[sample_col_name], verbose=True)

    plot_pca(rna.iloc[:, fea_start_id:], components=[1, 2], figsize=(8, 7),
             # color_vector = rna['CELL'].map(lambda x: src_name_map[x.split('.')[0].lower()]),
             color_vector=src_names_arr,
             scale=True, title='PCA (raw data)');
    plt.savefig(outpath/'pca_raw.png', dpi=dpi)

    # ---------------------------------------------------------
    #   Source scaling
    # ---------------------------------------------------------
    print('\nSource scaling ...')

    # Source scale
    rna_src = scale_rna(rna, fea_start_id=fea_start_id, per_source=True)
    print(rna_src.shape)

    plot_pca(rna_src.iloc[:, fea_start_id:], components = [1, 2], figsize=(8, 7),
             # color_vector=rna['CELL'].map(lambda x: src_name_map[x.split('.')[0].lower()]),
             color_vector=src_names_arr,
             scale=True, title='PCA (source-scaled data)');
    plt.savefig(outpath/'pca_src_scale.png', dpi=dpi)

    # Map gene expression to other cell lines
    df = map_rna_profiles(rna_src, cl_map=cl_mapping, sample_col_name=sample_col_name)
    df.to_csv(outpath/(datapath.name + '_src_scale'), index=False, sep='\t')

    # print(df.Sample.map(lambda x: x.split('.')[0].lower()).value_counts())
    # print(rna_copy.Sample.map(lambda x: x.split('.')[0].lower()).value_counts())
    # set(df.Sample.tolist()) - set(rna_copy.Sample.tolist())

    # ---------------------------------------------------------
    #   Combat scaling
    # ---------------------------------------------------------
    print('\nComBat scaling ...')

    # Data for combat
    meta = src_from_cell_col(rna[sample_col_name], verbose=False).to_frame()
    meta[sample_col_name] = rna[sample_col_name]

    # Combat scale
    rna_combat = combat_(rna, meta, sample_col_name=sample_col_name, batch_col_name='source')

    plot_pca(rna_combat.iloc[:, fea_start_id:], components=[1, 2], figsize=(8, 7),
             # color_vector = rna['CELL'].map(lambda x: src_name_map[x.split('.')[0].lower()]),
             color_vector=src_names_arr,
             scale=True, title='PCA (combat-scaled data)');
    plt.savefig(outpath/'pca_combat.png', dpi=dpi)

    # Map gene expression to other cell lines
    df = map_rna_profiles(rna_combat, cl_map=cl_mapping, sample_col_name=sample_col_name)
    df.to_csv(outpath/(datapath.name + '_combat'), index=False, sep='\t')

    print('Done.')


def main(args):
    args = parse_args(args)
    # args = vars(args)
    run(args)


if __name__ == '__main__':
    main(sys.argv[1:])


