#!/usr/bin/env python
# -*- coding: utf-8 -*- 

# @Author : Zongqin Wang 
# @File : FSQSSA.py
# @desc : A feature selection model based on quantum squirrel Search algorithm

import numpy as np
import math
import pandas as pd
import time
import os
from Cluster import Cluster


class FSQSSA:
    def __init__(self, data, data_label=None, k=None, max_iter=100, RESULT_PATH="./result"):
        '''
        Initialize the FSQSSA class and the required parameters
        :param data: scRNA-seq data
        :param data_label: cell label
        :param k: number of cluster
        :param max_iter: Maximum number of iterations
        :param RESULT_PATH: Result save path
        '''

        # Assigning values to variables
        self.max_iter = max_iter
        self.data = data
        self.cell = data.shape[0]
        self.gene = data.shape[1]
        self.k = k
        self.label = data_label

        # Clusterizer
        self.Cluster = Cluster(self.data, self.label, self.k)

        # FSQSSA requires various parameters
        self.Pdp = 0.1  # Probability of occurrence of natural enemies
        self.Gc = 2  # Sliding constants
        self.dg = lambda: np.random.uniform(0.3, 0.7)  # Random glide distance

        # Intermediate Result Variables
        self.sc_list = []  # SC score for each round
        self.fitness_list = []  # Fitness values for each round
        self.feature_list = []  # Number of features(genes) for each round
        self.subset_index = []  # Filtered feature index

        # Result save path
        self.RESULT_PATH = RESULT_PATH

    def random_location(self):
        '''
        The location of the squirrel's random escape
        :return: Returns a set of vectors representing the random position of a squirrel
        '''
        return np.random.uniform(0, 1, self.gene) * 2 * np.pi

    def init_population(self, squirrel_num=50):
        '''
        init population
        :param squirrel_num: Number of squirrels, default 50
        :return: Returns the initialized population
        '''
        pop = np.zeros([squirrel_num, self.gene], dtype=np.float64)
        for i in range(squirrel_num):
            pop[i] = self.random_location()
        return pop

    def run_task(self):
        '''
        FSQSSA run the task and the squirrel iterates max_iter times
        :return: Return gene subset
        '''
        population = self.init_population()
        for i in range(self.max_iter):
            population = self.task(population=population, epoch=i)

        self.cluster(population, self.max_iter)

        # save info
        info = {"fitness": self.fitness_list, "ari": self.sc_list, "feature": self.feature_list}

        t = time.time()
        if not os.path.exists(self.RESULT_PATH):
            os.makedirs(self.RESULT_PATH)

        subset = self.data[:, self.subset_index]

        info_path = self.RESULT_PATH + "/result.csv"
        subset_path = self.RESULT_PATH + "/subset.csv"
        subset_index_path = self.RESULT_PATH + "/subset_index.csv"

        pd.DataFrame(info).to_csv(info_path)
        pd.DataFrame(subset).to_csv(subset_path)
        pd.DataFrame(subset_index).to_csv(subset_index_path)

        # Return feature subset
        return subset

    def task(self, population, epoch):
        '''
        Single Iteration Tasks
        :param population population
        :param epoch epoch
        :return: Return a new population
        '''

        # Calculation Fitness
        indices = self.cluster(population, epoch)

        ht = indices[0]  # hickory nut tree
        at = indices[1:4]  # acorn nut trees
        nt = indices[4:]  # normal trees
        nt_1 = nt[[x for x in range(46) if x % 2 == 0]]  # Half flew to the acorn tree
        nt_2 = nt[[x for x in range(46) if x % 2 != 0]]  # Half flew to the hickory tree

        # Generate new locations
        new_pop = population.copy()

        # acorn nut trees (at) => hickory nut tree (ht)
        for index in at:
            if np.random.random() >= self.Pdp:
                new_pop[index] = new_pop[index] + self.dg() * self.Gc * (new_pop[ht] - new_pop[index])
            else:
                new_pop[index] = self.random_location()

        # normal trees (nt_1) => acorn nut trees (at)
        for index in nt_1:
            if np.random.random() >= self.Pdp:
                new_pop[index] = new_pop[index] + self.dg() * self.Gc * (new_pop[np.random.choice(at)] - new_pop[index])
            else:
                new_pop[index] = self.random_location()

        # normal trees (nt_2) => hickory nut tree (ht)
        for index in nt_2:
            if np.random.random() >= self.Pdp:
                new_pop[index] = new_pop[index] + self.dg() * self.Gc * (new_pop[ht] - new_pop[index])
            else:
                new_pop[index] = self.random_location()

        # Random relocation at the end of winter season
        s_min = 1e-5 / (365 ** ((epoch + 1) / (self.max_iter / 2.5)))
        sc = np.sqrt(np.sum((new_pop[at] - new_pop[ht]) ** 2))
        if sc < s_min:
            new_pop[nt_1] = 0 + self.levy_flight(size=population.shape[1]) * (2 * np.pi - 0)

        return new_pop

    def cluster(self, population, epoch):
        '''
        Print and save fitness, sc, feature, etc. for the current epoch
        :param population: Current population
        :param epoch: epoch
        :return: Return the index after sorting by fitness value
        '''
        fitness, ari, feature, subset_index_list = self.Cluster.fitness_kmeans(population)

        # Sort
        indices = np.argsort(-np.array(fitness))

        best_fitness = fitness[indices[0]]
        best_sc = ari[indices[0]]
        best_feature = feature[indices[0]]

        print("epoch {0} : fitness = {1} , sc = {2} , feature = {3}".format(
            epoch, best_fitness, best_sc, best_feature))

        self.fitness_list.append(best_fitness)
        self.feature_list.append(best_feature)
        self.sc_list.append(best_sc)
        self.subset_index = subset_index_list[indices[0]]

        return indices

    def levy_flight(self, alpha=0.01, beta=1.5, size=None):
        '''
        levy_flight Generate new location
        :return: new location
        '''
        sigma = (math.gamma(1 + beta) * np.sin(np.pi * beta / 2) / (
                math.gamma((1 + beta) / 2) * beta * 2 ** ((beta - 1) / 2))) ** (1 / beta)
        u = np.random.normal(0, sigma, size)
        v = np.random.normal(0, 1, size)
        sample = alpha * u / (np.abs(v) ** (1 / beta))
        return sample
