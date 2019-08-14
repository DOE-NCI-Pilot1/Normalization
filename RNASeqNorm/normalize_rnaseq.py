import os
import sys
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

file_path = os.path.dirname(os.path.relpath(__file__))
utils_path = os.path.abspath(os.path.join(file_path, '..', 'utils_py'))
sys.path.append(utils_path)

# todo: need to update pilot1_imports
from pilot1_imports import *

DATAPATH = '/Users/apartin/work/jdacs/Benchmarks/Data/Pilot1'
PDM_METADATA_FILENAME = 'combined_metadata_2018May.txt'
DEFAULT_DATATYPE = np.float16
SEED = 2018


def get_float_format(name):
    mapping = {}
    mapping['f16'] = np.float16
    mapping['f32'] = np.float32
    mapping['f64'] = np.float64

    mapped = mapping.get(name)
    if not mapped:
        raise Exception('No mapping found for "{}"'.format(name))
    return mapped


def init_params():
    parser = argparse.ArgumentParser(description='Normalize RNA-Seq data.')
    parser.add_argument('-out', '--outdir', dest='outdir',
                        default='.',
                        help='output dir to store the normalized rnaseq file')
    parser.add_argument('-f', '--in_fname', dest='in_fname',
                        default='combined_rnaseq_data_lincs1000',
                        help='rnaseq filename to normalize')
    parser.add_argument('-ff', '--float_format', dest='float_format',
                        default=argparse.SUPPRESS,
                        choices=['f16', 'f32', 'f64'],
                        help='float format of the output file')
    return parser.parse_args()


def run(args):
    print(args)
    dataset = args.in_fname
    outdir = args.outdir
    if ~hasattr(args, 'float_format'):
        float_format = DEFAULT_DATATYPE
    elif args.float_format in set(['f16', 'f32', 'f64']):
        get_float_format(args.float_format)


    # Load rnaseq
    # dataset = 'combined_rnaseq_data'
    # dataset = 'combined_rnaseq_data_lincs1000'
    df_rna = load_combined_rnaseq(dataset=os.path.join(DATAPATH, dataset), chunksize=2000, verbose=True)


    # Choose datasets to process
    # datasets_to_keep = ['ccle', 'ctrp', 'gdsc', 'ncipdm', 'gcsi', 'nci60']
    datasets_to_keep = ['ccle', 'gdc', 'nci60', 'ncipdm']
    use_metadata_file = True


    # Get the metadata
    if use_metadata_file:
        # Load metadata file
        meta = pd.read_csv(os.path.join(DATAPATH, PDM_METADATA_FILENAME), sep='\t')
        # meta = update_metadata_comb(meta)
        meta = update_metadata_comb_may2018(meta)
        meta = extract_specific_datasets(meta, datasets_to_keep=datasets_to_keep)
        df_rna, meta = update_df_and_meta(df_rna=df_rna, meta=meta, on='Sample')
    else:
        # If not using metadata file then use the source as meta
        meta = extract_specific_datasets(df_rna, datasets_to_keep=datasets_to_keep)
        meta['source'] = meta['Sample'].map(lambda x: x.split('.')[0].lower())  # add `source` col
        meta = meta[['Sample', 'source']]  # keep only `Sample` and `source` cols
        df_rna, meta = update_df_and_meta(df_rna=df_rna, meta=meta, on='Sample')


    # Load the cell line mappings
    cl_mapping = pd.read_csv(os.path.join(DATAPATH, 'cl_mapping'), sep='\t', header=None)
    cl_mapping.rename(columns={0: 'from_cell', 1: 'to_cell'}, inplace=True)


    # Plot PCA of raw data
    out_filename = dataset + '_raw'
    pca_obj, fig = plot_pca(df_rna.iloc[:, 1:],
                            color_vector=df_rna['Sample'].map(lambda x: x.split('.')[0].lower()),
                            components=[1, 2], to_scale=True, title='pca - raw data')
    fig.savefig('pca_' + out_filename + '.png', bbox_inches='tight')


    # ======  Source scaling  ======
    print('\nSource scaling...')
    df_rna_sc = scale_rnaseq(df=df_rna, per_source=True)
    out_filename = dataset + '_source_scale'

    pca_obj, fig = plot_pca(df_rna_sc.iloc[:, 1:],
                            color_vector=df_rna_sc['Sample'].map(lambda x: x.split('.')[0].lower()),
                            components=[1, 2], to_scale=True, title='pca - source scale')
    fig.savefig('pca_' + out_filename + '.png', bbox_inches='tight')

    df_rna_sc = copy_rna_profiles_to_cell_lines(df_rna_sc, cl_mapping)
    df_rna_sc = df_rna_sc.sort_values('Sample')

    df_rna_sc.to_csv(os.path.join(outdir, out_filename), sep='\t', float_format=float_format, index=False)



    # ======  Combat  ======
    print('\nCombat...')
    df_rna_be = ap_combat(df_rna, meta)
    out_filename = dataset + '_combat'

    pca_obj, fig = plot_pca(df_rna_be.iloc[:, 1:],
                            color_vector=df_rna_be['Sample'].map(lambda x: x.split('.')[0].lower()),
                            components=[1, 2], to_scale=True, title='pca - combat')
    fig.savefig('pca_' + out_filename + '.png', bbox_inches='tight')

    df_rna_be = copy_rna_profiles_to_cell_lines(df_rna_be, cl_mapping)
    df_rna_be = df_rna_be.sort_values('Sample')

    df_rna_be.to_csv(os.path.join(outdir, out_filename), sep='\t', float_format=float_format, index=False)


def main():
    args = init_params()
    run(args)


if __name__ == '__main__':
    main()