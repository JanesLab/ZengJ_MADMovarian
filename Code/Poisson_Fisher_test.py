#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  3 17:00:09 2021

@author: Catalina Alvarez
"""

import pandas as pd
from scipy.stats import chisquare, fisher_exact
import numpy as np
from scipy.stats import poisson


#%% Reading the raw data - 1 month
df= pd.read_csv('/data/Clonal_distribution/counts_events_5monthsWT.csv', sep=",")
path = '/figures_results/Clonal_distribution_results/5_monthsWT/'

#%% Calculing the expected and observed divisiones
# lmbd = 0.8165680473372781 #-> 1 month
# lmbd = 1.1284046692607004 #-> 2 months
lmbd = 0.7338709677419355 #-> 5 months

k = 0
#For individual models
models = np.unique(df['model'])                 
results = pd.DataFrame(columns = ['model', 'P_less_2', 'model_less_2', 'P_more_2', 'model_more_2'])
cutoff = 2
for j in models:
    model = j
    size = len(df[(df['model'] == j)])
    outliers = df[(df['divisions'] > cutoff) & (df['model'] == model)]
    observed_abovex = len(outliers)
    observed_belowx = size - observed_abovex
    expected_abovex = (1 - poisson.cdf(cutoff, lmbd))*size
    expected_belowx = size - expected_abovex     
    results.loc[k] = [model, expected_belowx, observed_belowx, expected_abovex, observed_abovex]
    k = k+1

model2 = 'All mice'
size2 = len(df)
outliers2 = df[(df['divisions'] > cutoff)]
observed_abovex2 = len(outliers2)
observed_belowx2 = size2 - observed_abovex2
expected_abovex2 = (1 - poisson.cdf(cutoff, lmbd))*size2
expected_belowx2 = size2 - expected_abovex2     
results.loc[k] = [model2, expected_belowx2, observed_belowx2, expected_abovex2, observed_abovex2]
  
results.to_csv("/figures_results/Poisson_results/5_monthsWT/expected_observed_results_5monthsWT.csv")

#%% Fisher exact test -> 5 months
table = pd.read_csv('/figures_results/Poisson_results/5_monthsWT/expected_observed_results_5monthsWT.csv', sep=",", index_col=0)

models = table['model']
results = pd.DataFrame(columns = ['model', 'oddsratio', 'p_value'])
j = 0
for i in range(len(models)): 
    contingence = np.array([[table.iloc[i,1], table.iloc[i,3]], [table.iloc[i,2], table.iloc[i,4]]])
    oddsratio, pvalue = fisher_exact(contingence)
    j = j+2
    results.loc[j] = [models[i], oddsratio, pvalue]
results.to_csv("/figures_results/5_monthsWT/Fisher_test_results_5months.csv")









