import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns

from pandas.api.types import is_string_dtype, is_numeric_dtype

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# ComBat imports
import patsy
from combat import *


fea_prfx_dct = {'rna': 'GE_', 'cnv': 'CNV_', 'snp': 'SNP_',
                'dsc': 'DD_', 'fng': 'FNG_'}


# def load_rna(datadir, logger=None, keep_cells_only=True, float_type=np.float32, impute=True):  
def load_rna(datapath, logger=None, gene_prfx=None, keep_cells_only=True, float_type=np.float32, impute=True):  
    """ Load RNA-Seq data. """
    # fname = 'combined_rnaseq_data_lincs1000'
    na_values = ['na', '-', '']
    if logger is None:
        print_func = print
    else:
        print_func = logger.info
        
    # if logger: logger.info('\nLoad RNA-Seq ... {datadir / fname}')
    print_func('\nLoad RNA-Seq ... {datadir / fname}')
    # rna = pd.read_csv(Path(datadir)/fname, sep='\t', low_memory=False, na_values=na_values, warn_bad_lines=True)
    rna = pd.read_csv(Path(datapath), sep='\t', low_memory=False, na_values=na_values, warn_bad_lines=True)
    # rna = rna.astype(dtype={c: float_type for c in rna.columns[1:]})  # Cast features
    # rna = rna.rename(columns={c: fea_prfx_dct['rna']+c for c in rna.columns[1:] if fea_prfx_dct['rna'] not in c}) # prefix rna gene names
    if gene_prfx is not None:
        rna = rna.rename(columns={c: gene_prfx+c for c in rna.columns[1:] if gene_prfx not in c}) # prefix rna gene names
    rna.rename(columns={'Sample': 'CELL'}, inplace=True)

    # if logger: logger.info(f'rna.shape {rna.shape}')
    print_func(f'rna.shape {rna.shape}')
    return rna


def scale_rna(df, fea_start_id, per_source=False):
    """ Scale df values and return updated df. """
    df = df.copy()

    if per_source:
        sources = df['CELL'].map(lambda x: x.split('.')[0].lower()).unique().tolist()
        for i, source in enumerate(sources):
            print('Scaling {}'.format(source))
            source_vec = df['CELL'].map(lambda x: x.split('.')[0].lower())
            source_idx_bool = source_vec.str.startswith(source)

            fea = df.loc[source_idx_bool, df.columns[fea_start_id:].values]
            fea_scaled = StandardScaler().fit_transform(fea)
            df.loc[source_idx_bool, fea_start_id:] = fea_scaled
    else:
        if is_numeric_dtype(df.iloc[:, fea_start_id:]):
            df.iloc[:, fea_start_id:] = StandardScaler().fit_transform(df.iloc[:, fea_start_id:])

    return df


def src_from_cell_col(cell_name_col, verbose=False):
    """ Takes column that contains sample names, extract the source name, and returns the arr of source names.
    Prints value_counts if verbose=True. """
    src_names_arr = cell_name_col.map(lambda x: x.split('.')[0].lower())
    src_names_arr.name = 'source'
    if verbose:
        print(src_names_arr.value_counts())
    return src_names_arr


def combat_(rna, meta, sample_col_name:str, batch_col_name:str):
    """
    sample_col_name : name of the column that contains the rna samples
    batch_col_name : name of the column that contains the batch values
    """
    rna_fea, pheno, _, _ = py_df_to_R_df(data=rna, meta=meta, sample_col_name=sample_col_name)
    # dat.columns.name = None
    # pheno.index.name = pheno.columns.name
    # pheno.columns.name = None

    mod = patsy.dmatrix("~1", data=pheno, return_type="dataframe")
    ebat = combat(data = rna_fea,
                  batch = pheno[batch_col_name],  # pheno['batch']
                  model = mod)

    df_rna_be = R_df_to_py_df(ebat, sample_col_name=sample_col_name)
    return df_rna_be


def R_df_to_py_df(data, sample_col_name):
    """ This is applied to the output of combat.py """
    return data.T.reset_index().rename(columns={'index': sample_col_name})


def py_df_to_R_df(data, meta, sample_col_name, filename=None, to_save=False, to_scale=False, var_thres=None):
    """ Convert python dataframe to R dataframe (transpose). """
    # Transpose df for processing in R
    data_r = data.set_index(sample_col_name, drop=True)
    data_r = data_r.T

    # This is required for R
    meta_r = meta.set_index(sample_col_name, drop=True)
    del meta_r.index.name
    meta_r.columns.name = sample_col_name

    return data_r, meta_r, data, meta



def plot_pca(df, components=[1, 2], figsize=(8, 6),
             color_vector=None, marker_vector=None,
             scale=False, grid=False, title=None, verbose=True):
    """
    Apply PCA to input df.
    Args:
        color_vector : each element corresponds to a row in df. The unique elements will be colored
            with a different color.
        marker_vector : each element corresponds to a row in df. The unique elements will be marked
            with a different marker.
    Returns:
        pca_obj : object of sklearn.decomposition.PCA()
        pca : pca matrix
        fig : PCA plot figure handle
    """
    if color_vector is not None:
        assert len(df) == len(color_vector), 'len(df) and len(color_vector) must be the same size.'
        n_colors = len(np.unique(color_vector))
        colors = iter(cm.rainbow(np.linspace(0, 1, n_colors)))

    if marker_vector is not None:
        assert len(df) == len(marker_vector), 'len(df) and len(marker_vector) shuold be the same size.'
        all_markers = ('o', 'v', 's', 'p', '^', '<', '>', '8', '*', 'h', 'H', 'D', 'd', 'P', 'X')
        markers = all_markers[:len(np.unique(marker_vector))]

    df = df.copy()

    # PCA
    if scale:
        X = StandardScaler().fit_transform(df.values)
    else:
        X = df.values

    n_components = max(components)
    pca_obj = PCA(n_components=n_components)
    pca = pca_obj.fit_transform(X)
    pc0 = components[0] - 1
    pc1 = components[1] - 1

    # Start plotting
    fig, ax = plt.subplots(figsize=figsize)

    if (color_vector is not None) and (marker_vector is not None):
        for i, marker in enumerate(np.unique(marker_vector)):
            for color in np.unique(color_vector):
                # print(i, 'marker:', marker, 'color:', color)
                idx = (marker_vector == marker) & (color_vector == color)
                ax.scatter(pca[idx, pc0], pca[idx, pc1], alpha=0.5,
                            marker=markers[i],
                            edgecolors='black',
                            color=next(colors),
                            label='{}, {}'.format(marker, color))

    elif (color_vector is not None):
        for color in np.unique(color_vector):
            idx = (color_vector == color)
            ax.scatter(pca[idx, pc0], pca[idx, pc1], alpha=0.5,
                        marker='o',
                        edgecolors='black',
                        color=next(colors),
                        label='{}'.format(color))

    elif (marker_vector is not None):
        for i, marker in enumerate(np.unique(marker_vector)):
            idx = (marker_vector == marker)
            ax.scatter(pca[idx, pc0], pca[idx, pc1], alpha=0.5,
                        marker=markers[i],
                        edgecolors='black',
                        color='blue',
                        label='{}'.format(marker))

    else:
        ax.scatter(pca[:, pc0], pca[:, pc1], alpha=0.5,
                   marker='s', edgecolors='black', color='blue')

    if title: ax.set_title(title)
    # ax.set_xlabel('PC' + str(components[0]))
    # ax.set_ylabel('PC' + str(components[1]))
    ax.set_xlabel('PC' + str(components[0]) + ' (explained var {:.3f})'.format(pca_obj.explained_variance_ratio_[pc0]))
    ax.set_ylabel('PC' + str(components[1]) + ' (explained var {:.3f})'.format(pca_obj.explained_variance_ratio_[pc1]))
    ax.legend(loc='lower left', bbox_to_anchor=(1.01, 0.0), ncol=1, borderaxespad=0, frameon=True)
    if grid:
        plt.grid(True)
    else:
        plt.grid(False)

    if verbose:
        print('Explained variance by PCA components [{}, {}]: [{:.5f}, {:.5f}]'.format(
            components[0], components[1],
            pca_obj.explained_variance_ratio_[pc0],
            pca_obj.explained_variance_ratio_[pc1]))

    plt.tight_layout()

    return pca_obj, pca, fig
