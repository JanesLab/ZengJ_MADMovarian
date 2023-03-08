#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 12:27:00 2021

@author: Catalina Alvarez
"""

import pandas as pd
import seaborn as sns
from scipy.stats import poisson, ks_2samp
import numpy as np
import matplotlib.pyplot as plt

#Reading the raw data
df= pd.read_csv('/data/Clonal_distribution/counts_events_5monthsWT.csv', sep=",")
path = '/figures_results/Clonal_distribution_results/5_monthsWT/'

#%% First method events -  KS test our data agains pmf
# models = ['#16816', '#16836', '#16935'] # 1 month
# models = ['#14548', '#17088', '#17089', '#17091'] # 2 months
# models = ['#14548', '#15655', '#15656'] # 5 months
models = np.unique(df['model'])

#Choose one model and evaluate the significance of the Poisson distribution
models = ['#14548']
for j in models:
    upper = [7]
    for z in upper:
        new_df = df[(df['divisions'] < z) & (df['model'] == j)]
        mydata = new_df['divisions']
        mlest = mydata.mean()
        lenx = len(mydata)
        k = np.arange(mydata.max()+1)
        k = k.astype(int)
        lenk = len(k)
        ipoisson = poisson.pmf(k, mlest)*lenx
        idealpoisson = ipoisson.astype(int)
        dpoisson = []
        text = []
        for i in range(lenk):
            dpoisson = dpoisson + [k[i]]*idealpoisson[i]
        stat, p_value = ks_2samp(mydata, dpoisson)       
        sns.set_style("white") 
        ax = sns.ecdfplot(mydata, label=j)
        ax = sns.ecdfplot(dpoisson , label='Poisson')
        plt.hist(text, label= 'Data mean = ' + str(round(mlest,2)) + '\n' + 'p_value = ' + str(round(p_value, 4)) + '\n' + 'KS = ' +  str(round(stat, 2)), color = 'white')
        ax.legend(loc="lower right")
        ax.set_title('Clonal events CDF ' + j)
        ax.set_xlabel('Cell division events', fontsize=12)
        ax.set_ylabel('CDF', fontsize=12)
        plt.savefig(path + j + '_ks_cdf_' + str(z) + '.png', dpi=300)
        plt.show()
     



































