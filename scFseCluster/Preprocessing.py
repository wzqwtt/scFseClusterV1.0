#!/usr/bin/env python
# -*- coding: utf-8 -*- 

# @Author : Zongqin Wang 
# @File : Preprocessing.py
# @desc : Read and preprocess scRNA-seq data

import pandas as pd
import numpy as np
import scanpy as sc
import os


def read_data(data_path, label_path=None, get_HVGs=True, k=None, save_path="./result"):
    '''
    Read and preprocess scRNA-seq data
    :param data_path: scRNA-data path
    :param label_path: label path
    :param get_HVGs: Acquisition of highly variable genes or not
    :param k: Number of clusters. It can be automatically evaluated by scFseCluster
    :param save_path: Result save path
    :return: data, label and k
    '''

    # read data
    data = pd.read_csv(data_path, index_col=0)
    print("data shape = {}".format(data.shape))

    # Selecting highly variable genes
    if get_HVGs == True:
        data, _ = Selecting_highly_variable_genes(data, 2000)
        # make dir
        if not os.path.exists(save_path):
            os.makedirs(save_path)
            print("make dir {}".format(save_path))
        # save the T2000 HVGs
        pd.DataFrame(data.T).to_csv(save_path + '/HVGs_expression_matrix.txt')
    else:
        data = np.array(data).T
    print("HVG expression matrix shape = {}".format(data.shape))

    # read label
    label = None
    if label_path is not None:
        label = pd.read_csv(label_path)
        label = np.array(label['cell_type'])
        k = len(np.unique(label))
        print("label is not none. label shape = {}, k = {}".format(label.shape, k))
    else:
        # Evaluate the number of K
        if k is None:
            print("Evaluate the number of K...")
            k = Evaluate_k(data)
        print("label is none. k = {}".format(k))

    return data, label, k


def Selecting_highly_variable_genes(X, highly_genes=2000):
    '''
    Selecting highly variable genes
    :param X: scRNA-seq data
    :param highly_genes: The number of highly variable genes. Default 2000
    :return: data, adata
    '''
    adata = sc.AnnData(X)
    adata.var_names_make_unique()
    # sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=highly_genes)
    adata = adata[:, adata.var['highly_variable']].copy()
    data = adata.X
    return data, adata


def Evaluate_k(data):
    '''
    Evaluate the number of K
    :param data: HVGs expression matrix
    :return: the number of cluster
    '''
    adata = sc.AnnData(data)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.louvain(adata)
    k = len(np.unique(adata.obs['louvain']))
    return k
