#!/usr/bin/env python
# -*- coding: utf-8 -*- 

# @Author : Zongqin Wang 
# @File : Cluster.py
# @desc : Evaluate the fitness of each feature subset

import numpy as np
from sklearn.cluster import KMeans
from sklearn import metrics


class Cluster:
    def __init__(self, data, data_label, k):
        '''
        Perform clustering tasks, mainly to obtain fitness values
        :param data: scRNA-seq data
        :param data_label: cell label
        :param k: the number of cluster
        '''

        self.data = data  # scRNA-seq data
        self.cell = data.shape[0]  # the number of cell
        self.gene = data.shape[1]  # the number of gene
        self.label = data_label  # cell label
        self.k = k

        self.random_state = 42  # Random number seeds
        self.w = 0.9  # Fitness value weights

        self.sc = []  # sc values for each squirrel
        self.fitness = []  # fitness values for each squirrel
        self.feature = []  # Number of features for each squirrel
        self.subset_index_list = []  # Index of features for each squirrel

    def get_feature_subset(self, squirrel):
        '''
        Find the subset of features corresponding to a quantum squirrel
        :param squirrel: A squirrel in the population
        :return: Returns the subset of features corresponding to this squirrel and the index corresponding to the features
        '''
        subset = []
        for i in range(squirrel.size):
            theta = squirrel[i]
            alpha = np.cos(theta) * np.cos(theta)  # Probability of 0
            beta = np.sin(theta) * np.sin(theta)  # Probability of 1
            if alpha < beta:
                subset.append(i)
        feature_subset = self.data[:, subset]
        # The following code preserves the genetic information
        # pd.DataFrame(feature_subset,columns=subset)
        return feature_subset, subset

    def fitness_without_label(self, subset, label, feature):
        '''
        Calculating fitness values
        :param subset: Current epoch gene subset
        :param label: kmeans output labels
        :param feature: the number of feature
        :return: return fitness and sc value
        '''
        f = (metrics.silhouette_score(subset, label, metric='euclidean') + 1) / 2
        if feature < 700: return -100, 0
        return self.w * f + (1 - self.w) * (1 - (feature / self.gene)), f

    def fitness_kmeans(self, population):
        '''
        Calculate the fitness of this population using the KMeans algorithm
        :param population: population
        :return: Returns the sc value, the number of features, and the fitness value; both are vectors
        '''

        self.sc.clear()
        self.fitness.clear()
        self.feature.clear()
        self.subset_index_list.clear()

        for unity in population:
            # Pick the subset of features corresponding to the individual with the feature index
            feature_subset, subset_index = self.get_feature_subset(unity)

            self.subset_index_list.append(subset_index)
            kmeans = KMeans(n_clusters=self.k, random_state=self.random_state).fit(feature_subset)  # Kmeans

            self.feature.append(feature_subset.shape[1])  # Add the number of feature subsets to the list
            fit, f = self.fitness_without_label(feature_subset, kmeans.labels_, feature_subset.shape[1])
            self.fitness.append(fit)
            self.sc.append(f)

        return self.fitness, self.sc, self.feature, self.subset_index_list
