#!/usr/bin/env python3

#Compute ASW over AnnData objects

import argparse
import numpy as np
import pandas as pd
import scanpy as sc 
import anndata


sc.settings.figdir = '.'

def silhouette(adata, group_key, metric='euclidean', embed='X_pca', scale=True):
    """
    wrapper for sklearn silhouette function values range from [-1, 1] with 1 being an ideal fit, 0 indicating overlapping clusters and -1 indicating misclassified cells
    """
    import sklearn.metrics as scm
    
    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f'{embed} not in obsm')
    asw = scm.silhouette_score(adata.obsm[embed], adata.obs[group_key])
    if scale:
        asw = (asw + 1)/2
    return asw

def compute_asw(adata):
    #apply function
    try:
        knn_graph = adata.uns['neighbors']
        print('BBKNN corrected object!')

    except KeyError:
        #Both: pre corrected logcounts and Scanorama counts enter through this way.
        #compute neighbors
        sc.tl.pca(adata, n_comps = args.n_pcs)
        if 'X_corrected' in adata.obsm.keys():
            rep_key = 'X_corrected'
        else:
            rep_key = 'X_pca'
        sc.pp.neighbors(adata, n_neighbors = args.n_neighbors, use_rep = rep_key, knn = True)
        
    sc.tl.umap(adata, n_components = args.n_pcs)
    sc.pl.umap(adata, color=[args.batch_key, args.celltype_key], legend_loc='on data', save=args.output_png, show=False)

    batch_asw = silhouette(adata, group_key = args.batch_key, embed = 'X_umap', scale = True)
    cell_asw = silhouette(adata, group_key = args.celltype_key, embed = 'X_umap', scale = True)
    print("ASW calculated!")
    
    results = pd.DataFrame(columns = ['Batch_ASW', 'Cell_type_ASW'])
    results.loc[0] = [batch_asw, cell_asw]

    save_file_to_csv(results)

def save_file_to_csv(results):
    results.to_csv(args.output_asw, header = True, index = False, sep = '\t')

def read_h5ad(file):
    
    #read input h5ad
    dataset =  sc.read(file)
    print("File read!")
    compute_asw(dataset)


# args
if __name__== "__main__":

    parser = argparse.ArgumentParser(description='Input/Output files')

    parser.add_argument("--input_object", 
                        dest='input_object',
                        type=str,
                        help ='Input h5ad object')
    parser.add_argument("--batch_key", 
                        dest='batch_key',
                        type=str,
                        default='Batch',
                        help ='Cell key defining Batch')
    parser.add_argument("--celltype_key", 
                        dest='celltype_key', 
                        type = str,
                        default= 'cell_type1', 
                        help ='Cell key defining cell type')
    parser.add_argument("--n_pcs", 
                        dest='n_pcs', 
                        type = int,
                        default= 25, 
                        help ='Number of PCs for PCA prior to graph correction')
    parser.add_argument("--n_neighbors", 
                        dest='n_neighbors', 
                        type= int,
                        default = 30,
                        help ='Number of nearest neighbours per batch to perform graph correction.')
    parser.add_argument('--output_asw', 
                        dest='output_asw',
                        type=str, 
                        help='Csv with ASW values')
    parser.add_argument('--output_png', 
                        dest='output_png',
                        type=str, 
                        help='Png of UMAPs')
    args = parser.parse_args()

    read_h5ad(args.input_object)


