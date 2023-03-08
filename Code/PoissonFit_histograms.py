#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 20:14:23 2021

@author: Catalina Alvarez
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import poisson
import numpy as np

#Reading the raw data
df= pd.read_csv('/data/Clonal_distribution/counts_events_5monthsWT.csv', sep=",")
path = '/figures_results/Clonal_distribution_results/5_monthsWT/'

#%% Poisson fit (pmf: probability mass function -> discrete probability distribution)
#https://stackoverflow.com/questions/37500406/how-to-fit-a-poisson-distribution-with-seaborn
#sns.distplot(x, kde = False, fit = stats.poisson)
# models = ['#16816', '#16836', '#16935'] # 1 month
# models = ['#14548', '#17088', '#17089', '#17091'] # 2 months
# models = ['#14548', '#15655', '#15656'] # 5 months
models = np.unique(df['model'])

#Choose one model and generate the histograms for the distribution taking a max. number of divisions (upper) as reference
models = ['#14548']
for j in models:
    upper = [14]
    for i in upper:
        new_df = df[(df['divisions'] < i) & (df['model'] == j)]
        x_final = new_df['divisions']
        mlest = x_final.mean()
        npoints = len(x_final)
        k = np.arange(x_final.max()+1)
        k = k.astype(int)
        # Random data
        # random = poisson.rvs(mlest, size=npoints)
        fig, ax1 = plt.subplots()
        ax = sns.countplot(x_final, order=k, alpha=0.7, color='r', label= j + ' ' + 'n < ' + str(i))
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=0)
        plt.plot(k, poisson.pmf(k, mlest)*len(x_final), 'bo', markersize=5, label='Poisson distribution')
        #ax = sns.countplot(random , alpha=0.7, color='g', label='Random data')
        plt.legend(loc="upper right")
        plt.title('Cell divisions in the Fallopian tube')
        plt.xlabel('Cell divisions', fontsize=12)
        plt.ylabel('Counts', fontsize=12)
        plt.setp(labels, rotation=0, fontsize=10)
        plt.savefig(path + j + '_Poissonfit_' + str(i) + '.png', dpi=300)
        plt.show()


