if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("M3C")
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/factoextra")
devtools::install_github("vqv/ggbiplot")
library("ggbiplot") 
library(factoextra)
library(M3C)
library(data.table)

#Reading our data (genes, samples TPMs) and creating a dataframe:
#Dataset for all samples (28: 12 Ray Samples, 4 FT, 12 TS) and 129 genes
data_all <- read.csv(file = "/data/DeSeq_GSEA_PCA/GSE135590_DEGs_TPMs_mapped.csv", row.names = 1)
data <- transpose(data_all)
colnames(data) <- rownames(data_all)
rownames(data) <- colnames(data_all)
data = data[7:34,]
matlog2 = log2(data+1)
mydata <- matlog2[1:12,]
gse <- matlog2[13:28,] 

#Running PCA
pca <- prcomp(gse) 
fviz_eig(pca)

#Standard plot for initial PCA
fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)     # Avoid text overlapping

#Costumized plots for initial PCA
fviz_pca_ind(pca, repel = TRUE, label = "none", title = "129 DEG in GSE dataset", invisible="quali", palette = c("steelblue2", "coral2"), habillage =c("GEO-FT", "GEO-FT", "GEO-FT", "GEO-FT","GEO-TS","GEO-TS", "GEO-TS", "GEO-TS", "GEO-TS", "GEO-TS", "GEO-TS", "GEO-TS", "GEO-TS", "GEO-TS", "GEO-TS", "GEO-TS"), labelsize = 10,
             pointsize = 2, ggtheme=theme(axis.text=element_text(size=10), axis.title=element_text(size=10), title=element_text(size=10)))

#Projection of our data into the GSE135590 PCA
scale <- scale(mydata, pca$center, pca$scale) %*% pca$rotation 

#Consolidating the PCA results
pca.plusproj.scale <- pca
pca.plusproj.scale$x <- rbind(pca.plusproj.scale$x, scale)
fviz_eig(pca.plusproj.scale)

#Standard plot for the final PCA
fviz_pca_ind(pca.plusproj.scale,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)     # Avoid text overlapping

#Costumized plots for final PCA
fviz_pca_ind(pca.plusproj.scale, repel = TRUE, label= FALSE, invisible="quali", geom.ind = "point", title = "129 DEG in GEO datasets", palette = c("steelblue2", "coral2", "green3", "violet"), habillage =c("GEO-FT", "GEO-FT", "GEO-FT", "GEO-FT","GEO-TS", "GEO-TS", "GEO-TS", "GEO-TS","GEO-TS", "GEO-TS", "GEO-TS", "GEO-TS", "GEO-TS", "GEO-TS", "GEO-TS", "GEO-TS", "LC", "SC", "LC", "SC", "LC", "SC", "LC", "SC", "LC",
                                                                                                                                                      "SC", "LC","SC"), labelsize = 2,
             pointsize = 2, ggtheme=theme(axis.text=element_text(size=10), axis.title=element_text(size=10), title=element_text(size=10)))

#Results summary:
summary <- summary(pca.plusproj.scale)
loadings <-pca.plusproj.scale$rotation
PCs <- pca.plusproj.scale$x
setwd("/figures_results/DeSeq_GSEA_PCA_results/")
write.csv2(loadings, file = "PCA_loadings.csv")

#Calculate centroids:
#For Fallopian tube samples:
GEO_FTx <- mean(pca.plusproj.scale$x[1:4,1])
GEO_FTy <- mean(pca.plusproj.scale$x[1:4,2])
#For Tumor samples:
GEO_TSx <- mean(pca.plusproj.scale$x[5:16,1])
GEO_TSy <- mean(pca.plusproj.scale$x[5:16,2])
#For small clones:
SCx <- mean(pca.plusproj.scale$x[c(18,20,22,24,26,28),1])
SCy <- mean(pca.plusproj.scale$x[c(18,20,22,24,26,28),2])
#For large clones:
LCx <- mean(pca.plusproj.scale$x[c(17,19,21,23,25,27),1])
LCy <- mean(pca.plusproj.scale$x[c(17,19,21,23,25,27),2])

#Distances for the large (LC) and small clones (SC) to the GEO centroids:
raysamples <- pca.plusproj.scale$x[17:28,1:2]
centroid_FT <- c(GEO_FTx, GEO_FTy)
centroid_TS <- c(GEO_TSx, GEO_TSy)
raysamples <- rbind(raysamples, centroid_FT)
raysamples <- rbind(raysamples, centroid_TS)
raysamples <- as.data.frame(raysamples)
distances <- dist(raysamples, method = "euclidean")
result_distance<-as.matrix(distances)
setwd("/figures_results/DeSeq_GSEA_PCA_results/")
write.csv2(result_distance, file = "Eucledian_distance_to_FT&TScentroids.csv")

#Distances for all the samples to all centroids:
all_sammples <- pca.plusproj.scale$x[1:28,1:2]
centroid_SC <- c(SCx, SCy)
centroid_LC <- c(LCx, LCy)
centroid_FT <- c(GEO_FTx, GEO_FTy)
centroid_TS <- c(GEO_TSx, GEO_TSy)
all_sammples <- rbind(all_sammples, centroid_SC)
all_sammples <- rbind(all_sammples, centroid_LC)
all_sammples <- rbind(all_sammples, centroid_FT)
all_sammples <- rbind(all_sammples, centroid_TS)
distances <- dist(all_sammples, method = "euclidean")
result_distance<-as.matrix(distances)
row.names(raysamples)[15] <- "centroid_SC" 
row.names(raysamples)[16] <- "centroid_LC" 
setwd("/figures_results/DeSeq_GSEA_PCA_results/")
write.csv2(result_distance, file = "Eucledian_distance_to_allcentroids.csv")
