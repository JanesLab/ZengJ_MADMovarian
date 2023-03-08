library("pheatmap")

#create matrix from dataframe of data (genes, logFC values)
# leading.edge.mat = as.matrix(leading.edge.df)

data <- read.csv(file = "/data/DeSeq_GSEA_PCA/DEGs_paired_cherryPick_FDR10%_TPMs.csv", header = TRUE, row.names = 1)
matlog2 = log2(data+1)
mat=data.matrix(data, rownames.force = NA)
colnames(mat) <- c("18996-S","18999-S","19151-S", "19152-S", "19153-S",	"19313-S", "18996-L","18999-L", "19151-L", "19152-L",	"19153-L", "19313-L")


#create function for z-scoring the data to normalize by gene
cal_z_score <- function(x) { 
  (x-mean(x)) / sd(x)
}

#z-score normalized matrix for heatmap
data_mat_norm = t(apply(mat,1,cal_z_score))
hist(data_mat_norm)

#plot heatmap using pheatmap package
pheatmap(data_mat_norm,
         color = rev(RColorBrewer::brewer.pal(11, "RdBu")),
         cluster_rows = TRUE,
         cluster_cols = FALSE, clustering_distance_rows = "euclidean",
         clustering_method = "ward.D",
         scale = "none",
         fontsize_col = 8, fontsize_row = 5,
         border_color = "NA",
         angle_col = 90,
         # xlab = "Samples", ylab = "DEG",
         main = "",
         legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, max(data_mat_norm)),
         # legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, 3.1),
         legend_labels = c("-3", "-2", "-1", "0", "1", "2", "-3","z-score \n\n"),
         # legend_labels = c("-2", "-1", "0", "1", "2", "Normalized TPMS\n"),
         # legend = TRUE,
         fontsize = 8
         ) 

#function to save heatmap
save_pheatmap_pdf <- function(x, filename, width=8, height=10) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(heatmap3,paste0(output_dir,"Heatmap_",DE_pair,".pdf"))