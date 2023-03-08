#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 29 23:23:28 2021

@author: Catalina Alvarez
"""

import pandas as pd
from scipy.stats import poisson, ks_2samp
import numpy as np
import matplotlib.pyplot as plt

#Reading the raw data
df= pd.read_csv('/data/Clonal_distribution/Clonal_distribution/counts_events_5monthsWT.csv', sep=",")
path = '/figures_results/Clonal_distribution_results/5_monthsWT/'

#%% Calculating parameters for the final distribution
lower = 0
upper = 7
df1_all= df[(df['divisions'] >= lower) & (df['divisions'] < upper)] 
df1 = df1_all['divisions']
lmbd = np.mean(df1)


#%% Generating the null distributions and evaluating significance
results = pd.DataFrame(columns = ['model', 'KS-stat', 'p_value'])
y = 0
# Individual models:
models = np.unique(df['model'])
#All data
# models = ['All_datasets']
for j in models:
    upper = [2]
    for z in upper:
        # Individual models:
        new_df = df[(df['divisions'] < z) & (df['model'] == j)]
        mydata = new_df['divisions']
        #All data:
        # new_df = df[(df['divisions'] < z)]
        # mydata = new_df['divisions']
        lenx = len(mydata)
        k = np.arange(mydata.max()+1)
        k = k.astype(int)
        lenk = len(k)
        # My data plot:
        #ax = sns.countplot(mydata, order=k, alpha=0.7, color='r', label='#10191')        
        
        ipoisson = poisson.pmf(k, lmbd)*lenx
        idealpoisson = ipoisson.astype(int)
        dpoisson = []
        text = []
        for i in range(lenk):
            dpoisson = dpoisson + [k[i]]*idealpoisson[i]
        ks_stat, p_value = ks_2samp(mydata, dpoisson)
        results.loc[y] = [j, ks_stat, p_value]
        y = y+1
        
        # The ideal Poisson data plot
        #ax = sns.countplot(dpoisson, order=k, alpha=0.7, color='r', label='Ideal Poisson')
        kss = []
        pss = []
        for i in range(1000):
            random = poisson.rvs(lmbd, size=lenx)
            ks_stat2, p_value2 = ks_2samp(random, dpoisson)
            kss = kss + [ks_stat2]
            pss = pss + [p_value2]        

        df_prob = pd.DataFrame({'col_one': kss})
        likelihood = df_prob[df_prob['col_one'] <= ks_stat]
        area_ks = np.trapz(likelihood['col_one'])
        area_total = np.trapz(kss)
        prob = 1 - ((area_ks)/(area_total))  
        plt.hist(kss, bins=20, alpha=0.45, label='n < ' + str(z), edgecolor = None, linewidth=0, density = True)
        plt.hist(text, label= 'p-value = ' + str(round(prob,5)) , color='white')
        plt.axvline(ks_stat, ls='--', color = 'red', label = 'Fallopian tube data')
        plt.title('KS Null distribution ' + j, fontsize=14)
        plt.legend(loc="upper right")
        plt.xlabel('KS statistics', fontsize=12)
        plt.ylabel('Frequency', fontsize=12)    
        plt.savefig(path + j + '_ks_null_fit_10cutoff' + str(z) + '.png', dpi=300)
        plt.show()



