#!/usr/bin/env python
# -*- coding: utf-8 -*-

# @Author : Zongqin Wang
# @File : scFseCluster.py
# @desc : Program execution entrance

import argparse
import pandas as pd
from Preprocessing import read_data
from FSQSSA import FSQSSA
from sklearn.cluster import KMeans
import os


def run(args):
    '''
    Run scFseCluster
    :param args: Parameters entered by the user
    '''

    # Parsing parameters
    data_path = args.data_path
    label_path = args.label_path
    get_HVGs = args.get_HVGs
    k = args.k
    max_iter = args.max_iter
    save_path = args.save_path
    do_plot = args.do_plot

    # read and preprocess data
    data, label, k = read_data(data_path, label_path, get_HVGs, k, save_path)

    # run FSQSSA
    print("=================== Run FSQSSA ===================")
    fsqssa = FSQSSA(data=data, max_iter=max_iter, data_label=label, k=k, RESULT_PATH=save_path)
    subset = fsqssa.run_task()

    # cluster
    kmeans = KMeans(n_clusters=k, random_state=42).fit(subset)  # Kmeans
    cluster_label = pd.DataFrame(kmeans.labels_, columns=["type"])
    cluster_label_path = save_path + "/cluster_label.csv"
    subset_path = save_path + "/subset.csv"

    # save file
    cluster_label.to_csv(cluster_label_path)

    # plot
    # Note: scFseCluster.py should be in the same y path as the R script plot_scFseCluster.R
    if do_plot:
        command = "Rscript plot_scFseCluster.R " + subset_path + " " + cluster_label_path
        os.system(command)


if __name__ == '__main__':
    # Program execution entrance
    parser = argparse.ArgumentParser(
        description="scFseCluster: a Feature selection enhanced clustering for single cell RNA-seq data.")

    parser.add_argument("-data", "--data_path", type=str, required=True,
                        help="scRNA-seq dataset path. You need to enter an absolute path or a relative path.")
    parser.add_argument("-label", "--label_path", type=str, default=None, required=False,
                        help="Cell labeling path for scRNA-seq dataset. You can choose whether you want to enter an "
                             "absolute path or a relative path.")
    parser.add_argument("-hvgs", "--get_HVGs", action="store_false", required=False,
                        help="When you do not use this parameter, it will help you filter HVGs. otherwise, it will "
                             "not help you filter.")
    parser.add_argument("-k", "--k", type=int, default=None, required=False,
                        help="You need to enter an int value specifying the number of clusters. scFseCluster can "
                             "automatically calculate the number of clusters, and by default, this parameter is not "
                             "required to be entered.")
    parser.add_argument("-iter", "--max_iter", type=int, default=100, required=False,
                        help="The number of FSQSSA iterations, default 100.")
    parser.add_argument("-plot", "--do_plot", action="store_true", required=False,
                        help="When you use this parameter, a scatter plot will be drawn.")
    parser.add_argument("-save", "--save_path", type=str, default="./result", required=False,
                        help="scFseCluster result output path, default output to the result directory under the "
                             "execution path.")

    args = parser.parse_args()
    print(args)

    # Run scFseCluster
    run(args)
