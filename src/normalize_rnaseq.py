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

from sklearn.impute import SimpleImputer, KNNImputer

# Utils
from utils import plot_pca, src_norm_rna, combat_norm

# Default settings
default_datapath = Path('./data/July2020/combined_rnaseq_data_lincs1000')
default_sources = ['ccle', 'nci60', 'ncipdm', 'gdc']
src_name_map = {'ccle': 'CCLE', 'nci60': 'NCI-60', 'ncipdm': 'NCI-PDM',
                'gdc': 'GDC', 'nova': 'Novartis'}
dpi = 100
float_format = np.float32

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
    # parser.add_argument('--include_gdc',
    #                     action='store_true',
    #                     help='Whether to include GDC (TCGA) data (default: False).')

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
        cl_map: df with columns 'from_sample' and 'to_sample' 
    Returns:
        df: contains the merged df which contains the original and the mapped samples
    """
    # Merge in order to copy rnaseq profiles
    mapped = cl_map.merge(rna, left_on='from_sample', right_on=sample_col_name, how='inner')
    
    # Drop and rename columns
    mapped = mapped.drop(columns=['from_sample', sample_col_name])
    mapped = mapped.rename(columns={'to_sample': sample_col_name})

    # Concat 'rna' (rna that contains original rna) and 'mapped' replicated rna
    df = pd.concat([rna, mapped], axis=0).sort_values(sample_col_name).reset_index(drop=True)
    return df


def run(args):
    import ipdb; ipdb.set_trace(context=5)
    fea_start_id = 1  # col index of the first gene feature
    datapath = Path(args.datapath)
    outpath = Path(args.outpath)
    fea_prfx = args.fea_prfx
    sample_col_name = 'Sample'
    os.makedirs(outpath, exist_ok=True)

    # ---------------------------------------------------------
    #   Load data
    # ---------------------------------------------------------
    print('\nLoad RNA-Seq ... {}'.format(datapath))
    rna = pd.read_csv(datapath, sep='\t', na_values=['na', '-', ''])
    if fea_prfx is not None:
        col_rename = {c: fea_prfx+c for c in rna.columns[1:] if fea_prfx not in c}
        rna = rna.rename(columns=col_rename)
    print(f'rna.shape {rna.shape}')

    # Load cell line mapping
    cl_mapping = pd.read_csv('data/cl_mapping', sep='\t', header=0,
                             names=['from_sample', 'to_sample'])

    # Keep or drop GDC (TCGA)
    # if args.include_gdc is False:
    #     aa = src_from_cell_col(rna[sample_col_name], verbose=False);
    #     idx = (aa != 'gdc')  # non-gdc indices
    #     rna = rna[idx].reset_index(drop=True)
    #     src_from_cell_col(rna[sample_col_name], verbose=False);

    # Drop duplicates in sample name
    # rna = rna.drop_duplicates(subset='Sample', keep=False)
    # rna = rna[ ~rna.Sample.isin(['GDC.BA2440R', 'GDC.BA2063R', 'GDC.BA2633R']) ].reset_index(drop=True)

    # ---------------------------------------------------------
    #   Extract sources that contain original RNA-Seq
    # ---------------------------------------------------------
    rna_copy = rna.copy()
    src_names_arr = src_from_cell_col(rna[sample_col_name], verbose=True)
    rna = rna_copy.loc[ src_names_arr.isin(args.src) ].reset_index(drop=True)
    print(f'rna.shape {rna.shape}')

    # ---------------------------------------------------------
    #   Impute missing values
    # ---------------------------------------------------------
    print('\nImpute missing values.')

    # Dump the count of NANs in each column
    aa = rna.isna().sum(axis=0).reset_index()
    aa = aa.rename(columns={'index': 'col', 0: 'count'})
    aa = aa.sort_values('count', ascending=False).reset_index(drop=True)
    aa.to_csv(outpath/'nan_count_per_col.csv', index=False)

    # Impute NaN values
    if sum(rna.isna().sum() > 0):
        print('Columns with NaNs: {}'.format( sum(rna.iloc[:, fea_start_id:].isna().sum() > 0) ))
        imputer = SimpleImputer(missing_values=np.nan, strategy='mean')
        # imputer = KNNImputer(missing_values=np.nan, n_neighbors=5,
        #                      weights='uniform', metric='nan_euclidean',
        #                      add_indicator=False)
        rna.iloc[:, fea_start_id:] = imputer.fit_transform(rna.iloc[:, fea_start_id:].values)
        print('Columns with NaNs: {}'.format( sum(rna.iloc[:, fea_start_id:].isna().sum() > 0) ))

        # Impute per source
        # TODO. There is a problem: 
        # src_names_arr = src_from_cell_col(rna[sample_col_name], verbose=True);
        # for s in src_names_arr.unique():
        #     idx = (src_names_arr == s).values
        #     imputer = SimpleImputer(missing_values=np.nan, strategy='mean')
        #     rna.iloc[idx, 1:] = imputer.fit_transform(rna.iloc[idx, 1:].values)

    # ---------------------------------------------------------
    #   Raw data
    # ---------------------------------------------------------
    print('\nRaw Data ...')
    src_names_arr = src_from_cell_col(rna[sample_col_name], verbose=True)

    plot_pca(rna.iloc[:, fea_start_id:], components=[1, 2], figsize=(8, 7),
             # color_vector=src_names_arr,
             color_vector=np.array([src_name_map[s] for s in src_names_arr]),
             scale=True, title='PCA (raw data)');
    plt.savefig(outpath/'pca_raw.png', dpi=dpi)

    # ---------------------------------------------------------
    #   Source Normalization
    # ---------------------------------------------------------
    print('\nSource Normalization ...')

    # Source scale
    rna_src = src_norm_rna(rna, fea_start_id=fea_start_id)
    rna_src = rna_src.astype({c: float_format for c in rna_src.columns[fea_start_id:]})
    print(rna_src.shape)

    plot_pca(rna_src.iloc[:, fea_start_id:], components = [1, 2], figsize=(8, 7),
             # color_vector=src_names_arr,
             color_vector=np.array([src_name_map[s] for s in src_names_arr]),
             scale=True, title='PCA (source-scaled data)');
    plt.savefig(outpath/'pca_src_scale.png', dpi=dpi)

    # Map gene expression to other cell lines
    df = map_rna_profiles(rna_src, cl_map=cl_mapping, sample_col_name=sample_col_name)
    df.to_csv(outpath/(datapath.name + '_src_scale'), index=False, sep='\t')

    # Is there any difference in samples?
    # print(df.Sample.map(lambda x: x.split('.')[0].lower()).value_counts())
    # print(rna_copy.Sample.map(lambda x: x.split('.')[0].lower()).value_counts())
    # set(df.Sample.tolist()) - set(rna_copy.Sample.tolist())

    # ---------------------------------------------------------
    #   Combat Normalization
    # ---------------------------------------------------------
    print('\nComBat Normalization ...')

    # Data for combat
    meta = src_from_cell_col(rna[sample_col_name], verbose=False).to_frame()
    meta[sample_col_name] = rna[sample_col_name]

    # Sample names must not have duplicates
    assert rna[sample_col_name].nunique() == rna.shape[0], \
        "Sample names must not have duplicates. Otherwise, ComBat results in error."

    # Combat scale
    rna_combat = combat_norm(rna, meta, sample_col_name=sample_col_name,
                             batch_col_name='source')
    rna_combat = rna_combat.astype({c: float_format for c in rna_combat.columns[fea_start_id:]})

    plot_pca(rna_combat.iloc[:, fea_start_id:], components=[1, 2], figsize=(8, 7),
             # color_vector=src_names_arr,
             color_vector=np.array([src_name_map[s] for s in src_names_arr]),
             scale=True, title='PCA (combat-scaled data)');
    plt.savefig(outpath/'pca_combat.png', dpi=dpi)

    # Map gene expression to other cell lines
    df = map_rna_profiles(rna_combat, cl_map=cl_mapping, sample_col_name=sample_col_name)
    df.to_csv(outpath/(datapath.name + '_combat'), index=False, sep='\t')

    # ---------------------------------------------------------
    print('Done.')


def main(args):
    args = parse_args(args)
    run(args)


if __name__ == '__main__':
    main(sys.argv[1:])
