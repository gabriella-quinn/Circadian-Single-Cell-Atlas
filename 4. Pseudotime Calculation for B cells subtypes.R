###############################################################
###############################################################
#      R script corresponding to Maples et al.                #
# A single-cell atlas of intestinal immune cells across the   # 
#   day-night cycle reveals dynamic populations               #
###############################################################
###############################################################

##All packages used in this script can found cited in the Key Resource Table attached to the manuscript##

##########################FOURTH SCRIPT OF FOUR#############################################
######################### Pseudotime Analysis of B cells####################################


####ISOLATE JUST B CELLS OUT OF THE LARGER SEURAT OBJECT 

rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(reticulate)
options(future.globals.maxSize = 300000 * 1024^2)

#loading in annotated single cell object 
#see "1. Integration and Annotation of Circacdian Small Intestine Atlas.R" for parameters used. 
load('IntegratedALLCELLS.RData')

##use the representative cells as calculated in "1. Integration and Annotation of Circacdian Small Intestine Atlas.R"
object <- subset(object, subset = sketch_snn_res.0.1 %in% c(0,1:10))

#isolate B cells for the lineages
bs <- subset(object, subset = annotations %in% keep)
table(bs$annotations)
DefaultAssay(bs) <- "RNA"
bs <- JoinLayers(bs)
annotations.for.later <- bs$annotations

##remake Seurat object to remove sketch assay and re-scale
bs <- CreateSeuratObject(counts = bs[['RNA']]$counts)
bs$annotations <- annotations.for.later
bs[["percent.mt"]] <- PercentageFeatureSet(bs, pattern = "^mt-")
bs <- subset(bs, subset = percent.mt < 2.5)
bs <- NormalizeData(bs, normalization.method = "LogNormalize", scale.factor = 10000)
bs <- FindVariableFeatures(bs, verbose = F)
bs <- ScaleData(bs, verbose = F)
bs <- RunPCA(bs, verbose = F, reduction.name = 'bs.pca')
bs <- RunUMAP(bs, reduction = "bs.pca", dims = 1:10, return.model = T, verbose = F)

bs <- FindNeighbors(bs, dims = 1:2, reduction = 'umap')
bs <- FindClusters(bs, resolution = 0.1)






####RUN MONOCLE3 TO CALCULATE PSEUDOTIME VALUE FOR EACH CELL

library(monocle3)
library(tidyverse)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)

cds <- as.cell_data_set(bs)
cds <- estimate_size_factors(cds)
head(colData(cds))
fData(cds)
cds@clusters@listData[["UMAP"]][["partitions"]] <- bs$seurat_clusters

#set active ident 
bs <- SetIdent(bs, value = "seurat_clusters") #res0.1
list.cluster <- bs@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- bs@reductions$umap@cell.embeddings

cds <- learn_graph(cds, use_partition = F)
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[,clusters(cds) == 1])) #Transitional B cell cluster

##identify genes that predict pseudotime value
deg <- graph_test(cds, neighbor_graph = "principal_graph")
deg %>% arrange(q_value) %>% filter(status == "OK") %>% head()

