# DESeq2 analysis for the Fallopian tube RNA data
library(BiocManager)
library("DESeq2")
library("tximport")
library("readr")
library("tximportData")
library(edgeR)
library(data.table)
library("apeglm")
library("ashr")
library("pheatmap")
library("vsn")
library(grid)
library(gridtext)
library("gplots")

# RSEM read counts were converted to transcripts per million (TPM) and used as input for this analysis
# Read the TPM data (ensgene, entrez, symbol, chr, ...) for all your samples and the samples information table
tpms <- read.csv(file = "/data/DeSeq_GSEA_PCA/rsem_counts2TPMs.csv", row.names = 3, header = TRUE)
tpms <- tpms[,9:20] # Select just the columns associated with the TPM per sample
tpms <- round(tpms) # Round the TPM values 

# Create and use a table of sample information
colnames(tpms) <- c("X18996_L", "X18996_S", "X18999_L", "X18999_S", "X19151_L", "X19151_S", "X19152_L", "X19152_S","X19153_L", "X19153_S", "X19313_L", "X19313_S") #Rename the columns accordingly with your samples (1st factor)
samples <- as.list(c(colnames(tpms))) #Generate a list with the samples names
conditions <- as.list(c("L", "S", "L", "S", "L", "S", "L", "S", "L", "S", "L", "S")) #Generate a list with the samples conditions (1st factor)
type <- as.list(c("sample1", "sample1", "sample2", "sample2", "sample3", "sample3", "sample4", "sample4", "sample5", "sample5", "sample6", "sample6")) #Generate a list with the samples source (types/2nd factor)
coldata <- cbind(samples, conditions, type) #Combine the factors as a list
coldata2 = as.data.frame(coldata) #Generate a dataframe with your factors
fwrite(coldata2,"coldataf_paired_final.csv") #Export the factors matrix
coldataf <- read.csv(file = "/data/DeSeq_GSEA_PCA/coldataf.csv", sep=",",row.names= 1, header = TRUE) #Read the factors matrix from a file

# Construct a DESeqDataSet:
dds <- DESeqDataSetFromMatrix(countData = tpms,
                              colData = coldataf,
                              design = ~ conditions)

# Multi-factor designs
ddsMF <- dds
levels(ddsMF$conditions) #Check levels in conditions
ddsMF$type <- factor(ddsMF$type, levels = c("sample1", "sample2", "sample3", "sample4","sample5", "sample6"))
levels(ddsMF$type) #Check levels in types
design(ddsMF) <- formula(~ type + conditions)

#Run the DEGs analysis and read the results
ddsMF <- DESeq(ddsMF)
resMF <- results(ddsMF) 
summary(resMF)

#The significance cutoff used for optimizing the independent filtering (by default 0.1). If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.
resOrdered <- resMF[order(resMF$pvalue),] #Take the list of results and order it by p-value
sum(resMF$padj < 0.1, na.rm=TRUE) #Identify how many genes have a padj-values <0.1
head(resMF) #Head of the total file
head(resOrdered) #Head of the padj.value ordered file
summary(resOrdered)

# Export data
write.csv(as.data.frame(resOrdered), 
          file="DEGs_paired_132genes_FDR10%.csv")

# Select different significance level (alpha) for adjusted FDR p-adj calculation
resMF2 <- results(ddsMF, alpha=0.2)
resOrdered2 <- resMF2[order(resMF2$pvalue),]
sum(resMF2$padj < 0.2, na.rm=TRUE)
summary(resMF2)

# Export data
write.csv(as.data.frame(resOrdered2), 
          file="DEGs_paired_132genes_FDR20%.csv")


