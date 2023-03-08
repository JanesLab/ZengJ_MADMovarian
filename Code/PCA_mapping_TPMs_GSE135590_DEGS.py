#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 11:28:09 2023

@author: catalinaalvarez
"""

import pandas as pd

# Import the GSE135590 file and get the list of all the genes in that dataset
database = pd.read_excel('/data/DeSeq_GSEA_PCA/GSE135590_TPMs.xlsx')
genes_db = database["symbol"].tolist()

# Import the file with the DEGs identified in our study and get the list of the genes. Before loading add the column names for the samples you want to map
degs = pd.read_excel('/data/DeSeq_GSEA_PCA/DEGs_paired_132genes_FDR10%.xlsx')
genes_degs = degs["genes"].tolist()

# Create the output file
writer = pd.ExcelWriter('GSE135590_DEGs_TPMs_mapped.xlsx', engine='xlsxwriter')

#Find and append the TPMs data from the GSE135590 file in our DEGs file
a=0
for i in range(len(genes_degs)):
    for j in range(len(genes_db)):
        if genes_degs[i]==genes_db[j]:
            degs['FT1'].iloc[i] = database['FT1'][j]
            degs['FT2'].iloc[i] = database['FT2'][j]
            degs['FT3'].iloc[i] = database['FT3'][j]
            degs['FT4'].iloc[i] = database['FT4'][j]
            degs['13604T'].iloc[i] = database['13604T'][j]
            degs['2611T'].iloc[i] = database['2611T'][j]
            degs['14576T'].iloc[i] = database['14576T'][j]
            degs['2942T'].iloc[i] = database['2942T'][j]
            degs['14109T'].iloc[i] = database['14109T'][j]
            degs['14104T'].iloc[i] = database['14104T'][j]
            degs['13179T'].iloc[i] = database['13179T'][j]
            degs['14397T'].iloc[i] = database['14397T'][j]
            degs['14474T'].iloc[i] = database['14474T'][j]
            degs['13085T'].iloc[i] = database['13085T'][j]
            degs['13576T'].iloc[i] = database['13576T'][j]
            degs['2163T'].iloc[i] = database['2163T'][j]

# Write the results to the output file
degs.to_excel(writer)
writer.save()
 