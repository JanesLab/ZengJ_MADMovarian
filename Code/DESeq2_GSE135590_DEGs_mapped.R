# DESeq2 analysis for GSE135590 RNA data mapped to our DEGs
library(BiocManager)
library("readr")
library(data.table)
library("DESeq2")

# RSEM read TPMs were downloaded from GEO database, mapped to the DEGs identified in our study (PCA_mapping_TPMs_GSE135590_DEGS.py) and used as input for this analysis
# Read the TPM data for all your samples and 129 DEGs
tpms <- read.csv(file = "/data/DeSeq_GSEA_PCA/GSE135590_DEGs_TPMs_mapped.csv", row.names = 1, header = TRUE)
tpms <- tpms[,23:34] #Select just the columns associated with the TPM per sample
tpms <- round(tpms) #Round the TPM values 

# Create and use a table of sample information
samples <- colnames(tpms) # Extract columns names  
samples <- as.list(samples) #Generate a list with the samples names (1st factor)
conditions <- as.list(c("co-clustered tumors", "co-clustered tumors","co-clustered tumors","co-clustered tumors","un-clustered tumors","co-clustered tumors","co-clustered tumors", "un-clustered tumors","un-clustered tumors", "un-clustered tumors","un-clustered tumors", "un-clustered tumors")) #Generate a list with the samples conditions (1st factor)
coldata <- cbind(samples, conditions) #Combine the factors as a list
coldata2 = as.data.frame(coldata) #Generate a dataframe with your factors
fwrite(coldata2,"GSE135590_coldataf.csv") #Export the factors matrix
coldataf <- read.csv(file = "/data/DeSeq_GSEA_PCA/GSE135590_coldataf.csv", sep=",",row.names= 1, header = TRUE) #Read the factors matrix from a file

# Construct a DESeqDataSet:
dds <- DESeqDataSetFromMatrix(countData = tpms,
                              colData = coldataf,
                              design = ~ conditions)

#Run the DEGs analysis and read the results
dds <- DESeq(dds)
res <- results(dds)
summary(res)
head(res) #Head of the total file

#The significance cutoff used for optimizing the independent filtering (by default 0.1). If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.
resOrdered <- res[order(res$pvalue),] #take the list of results and order it by p-value
sum(res$padj < 0.1, na.rm=TRUE) #identify how many genes have a padj-values <0.1
head(resOrdered) #Head of the padj.value ordered file
summary(resOrdered)

# Export data
write.csv(as.data.frame(resOrdered), 
          file="DEGs_GSE135590_16genes_FDR10%.csv")

# Select different significance level (alpha) for adjusted FDR p-adj calculation
res2 <- results(dds, alpha=0.2)
resOrdered2 <- res2[order(res2$pvalue),]
sum(res2$padj < 0.2, na.rm=TRUE)
summary(res2)

# Export data
write.csv(as.data.frame(resOrdered2), 
          file="DEGs_GSE135590_16genes_FDR20%.csv")


