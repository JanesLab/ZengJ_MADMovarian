#GSEA
#BiocManager::install("fgsea", force = TRUE)
library(data.table)
library(fgsea)

#The list of DEGs identified with DeSeq2 was used to map mouse genes to human orthologous (biomartRt.R). 
#Then, this list of genes was use to consolidate (GSEA_mapping_hGenes.py) the input excel file used here [gene names, padj, and FC[].

# data = data frame containing gene names, padj and FC
input <- read.csv(file = "/data/DeSeq_GSEA_PCA/gsea_input_paired_FDR10%.csv", header = TRUE)
# input <- read.csv(file = "/data/DeSeq_GSEA_PCA/gsea_input_GSE135590_allgenes.csv", header = TRUE)

data <- input[,c(4,5,6)] # Select just the columns associated with the human-names, padj and FC

#Load the msigb database
msigdb.url = 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.4/'

#GMT file options, choose one:
GO_file= paste0(msigdb.url, "msigdb.v7.4.symbols.gmt")
# GO_file= paste0(msigdb.url, "h.all.v7.4.symbols.gmt")

myGO = gmtPathways(GO_file)

# gene.list is ordered (by log2FC and padj) list of log2FC and gene names
gene.list <- setNames(data[order(data$log2FoldChange*-log(data$padj)),]$log2FoldChange,  data[order(data$log2FoldChange*-log(data$padj)),]$hsapiens_homolog_associated_gene_name)
gene.list = gene.list[!duplicated(names(gene.list))]
gene.list = gene.list[!is.na(gene.list)]

fgRes = fgsea(pathways = myGO,
                            stats = gene.list,
   
                                       minSize=5,
   
                                       maxSize=600)
fwrite(fgRes, file="gseaRes_msigdb_paired_FDR10%.txt", sep="\t", sep2=c("", " ", ""))
fwrite(fgRes, file="gseaRes_hallmarks_paired_FDR10%.txt", sep="\t", sep2=c("", " ", ""))

#fwrite(fgRes, file="gseaRes_msigdb_all_GSE135590.txt", sep="\t", sep2=c("", " ", ""))
#fwrite(fgRes, file="gseaRes_hallmarks_all_GSE135590.txt", sep="\t", sep2=c("", " ", ""))

#Plot
topPathwaysUp <- fgRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(myGO[topPathways], gene.list, fgRes, 
              gseaParam=0.5)

