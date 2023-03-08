#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd

# Import the DeSeq2 output file and get the list of all the genes in that dataset
degs = pd.read_excel('/data/DeSeq_GSEA_PCA/DEGs_paired_132genes_FDR10%.xlsx')
#degs = pd.read_excel('/data/DeSeq_GSEA_PCA/DEGs_GSE135590_1095genes_FDR10%.xlsx')
genes_degs = degs["genes"].tolist() 

# Import the human homologs list (biomartRt.R output file) and get the list of all the genes 
hgenes = pd.read_excel('/data/DeSeq_GSEA_PCA/filtered_orth_paired_FDR10%.xlsx')
#hgenes = pd.read_excel('/data/DeSeq_GSEA_PCA/filtered_orth_all_GSE135590.xlsx')
genes_hgenes = hgenes["external_gene_name"].tolist()

# Create the output file
writer = pd.ExcelWriter('gsea_input_paired_FDR10%.xlsx', engine='xlsxwriter')
#writer = pd.ExcelWriter('gsea_input_GSE135590_allgenes.xlsx', engine='xlsxwriter')

#Find and append the DeSeq2 data in our DEGs file
a=0
for i in range(len(genes_hgenes)):
    for j in range(len(genes_degs)):
        if genes_hgenes[i]==genes_degs[j]:
            hgenes['padj'].iloc[i] = degs['padj'][j]
            hgenes['log2FoldChange'].iloc[i] = degs['log2FoldChange'][j]

# Write the results to the output file            
hgenes.to_excel(writer)
writer.save()
