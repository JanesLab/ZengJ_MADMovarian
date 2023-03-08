# This is a handy function to get gene human homologs for mice genes
# Arg gene_vec: A vector of mice genes
# Arg mouse_biomart: mouse biomart db from which gene attributes will be extracted
# Value A list where the first item is a data.frame object with
# mice gene names, corresponding human homologs, chromosome and type of orthology.
# Second obj in the list is a vector with the names of genes that have no human homolog and
# thus could not be converted.
# BiocManager::install("biomaRt")
library(magrittr)
library(dplyr)
library("biomaRt")

browseVignettes("biomaRt")
#Read the file with the mouse gene list from our analysis
genelist <- read.csv(file = "/data/DeSeq_GSEA_PCA/DEGs_paired_132genes_FDR10%.csv", header = TRUE)
# genelist <- read.csv(file = "/data/DeSeq_GSEA_PCA/DEGs_GSE135590_1095genes_FDR10%.csv", header = TRUE)
gene_vec <- genelist[,1] #Select just the columns associated with the gene names

#Use the database to make the conversion
mouse_biomart = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

convert_mouse_to_human = function(gene_vec, mouse_biomart) {
  
  # Get attributes of interest from mouse biomart object
  mouse_attributes = c("external_gene_name",
                         "hsapiens_homolog_associated_gene_name",
                         "hsapiens_homolog_chromosome",
                         "hsapiens_homolog_orthology_type",
                         "hsapiens_homolog_perc_id")
    
  # Get human mouse orthologs
  # This is getting all of the orthologs. First filter for genes that have mouse
  # orthologs
  mouse_human_orthologs = getBM(attributes = mouse_attributes,
                                  filters = "with_hsapiens_homolog",
                                  values = TRUE,
                                  mart = mouse_biomart,
                                  uniqueRows = TRUE)
    
  # Filter for genes that have one-to-one human-mouse orthologs
  filtered_orth = mouse_human_orthologs %>%
     filter(hsapiens_homolog_orthology_type == "ortholog_one2one")
    
  # Select user genes that are one-to-one orthologs
  ol_idx = which(gene_vec %in% filtered_orth$external_gene_name)
    user_one_to_one = gene_vec[ol_idx]
    
  # Would be interesting know how many couldn't get converted
  unconverted_genes = setdiff(gene_vec, user_one_to_one)
    
    # Filter one-to-one ortholog list
  filtered_orth = filtered_orth %>%
      filter(external_gene_name %in% user_one_to_one)
    
    return(list(filtered_orth, unconverted_genes))
}

results = convert_mouse_to_human(gene_vec, mouse_biomart)
write.csv(as.data.frame(results[1]), 
          file="filtered_orth_paired_FDR10%.csv")
          # file="filtered_orth_all_GSE135590.csv")

write.csv(as.data.frame(results[2]), 
          file="unconverted_genes_paired_FDR10%.csv")
          # file="unconverted_genes_all_GSE135590.csv")
