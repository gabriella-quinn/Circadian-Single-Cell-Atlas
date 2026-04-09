###############################################################
###############################################################
#      R script corresponding to Maples et al.                #
# A single-cell atlas of intestinal immune cells across the   # 
#   day-night cycle reveals dynamic populations               #
###############################################################
###############################################################

##All packages used in this script can found cited in the Key Resource Table attached to the manuscript##

#######FIRST SCRIPT OF FOUR###########################################
#######Integrating 4 time points into a single Seurat object##########


rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(loupeR)
# set this option when analyzing large datasets
options(future.globals.maxSize = 3000000 * 1024^2)
set.seed(1234)


# ##loading in experimental data from circadian time course
zt0_2 <- Read10X('/ZT0_2/sample_filtered_feature_bc_matrix')
zt0_2 <- CreateSeuratObject(zt0_2, min.cells = 100)
zt0_2[["percent.mt"]] <- PercentageFeatureSet(zt0_2, pattern = "^mt-")
zt0_2 <- subset(zt0_2, subset = nFeature_RNA > 100 & percent.mt < 10)
zt0_2$sample <- 'ZT0'

zt6 <- Read10X('/ZT6/sample_filtered_feature_bc_matrix')
zt6 <- CreateSeuratObject(zt6, min.cells = 100)
zt6[["percent.mt"]] <- PercentageFeatureSet(zt6, pattern = "^mt-")
zt6 <- subset(zt6, subset = nFeature_RNA > 100 & percent.mt < 10)
zt6$sample <- 'ZT6'

zt12 <- Read10X('/ZT12/sample_filtered_feature_bc_matrix')
zt12 <- CreateSeuratObject(zt12, min.cells = 100)
zt12[["percent.mt"]] <- PercentageFeatureSet(zt12, pattern = "^mt-")
zt12 <- subset(zt12, subset = nFeature_RNA > 100 & percent.mt < 10)
zt12$sample <- 'ZT12'

zt18 <- Read10X('/ZT18/sample_filtered_feature_bc_matrix')
zt18 <- CreateSeuratObject(zt18, min.cells = 100)
zt18[["percent.mt"]] <- PercentageFeatureSet(zt18, pattern = "^mt-")
zt18 <- subset(zt18, subset = nFeature_RNA > 100 & percent.mt < 10)
zt18$sample <- 'ZT18'

object <- merge(zt0_2, y = c(zt6,zt12,zt18), add.cell.ids = c('zt0','zt6','zt12','zt18'))
object[["RNA"]] <- JoinLayers(object[["RNA"]])
object <- NormalizeData(object)
# split assay by sample
object[["RNA"]] <- split(object[["RNA"]], f = object$sample)
object <- FindVariableFeatures(object, verbose = FALSE)
#take a 5th of the cells to do the estimations on 
object <- SketchData(object = object, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")

DefaultAssay(object) <- "sketch"
object <- FindVariableFeatures(object, verbose = F)
object <- ScaleData(object, verbose = F)
object <- RunPCA(object, verbose = F)
# integrate the datasets
object <- IntegrateLayers(object, method = RPCAIntegration, orig = "pca", new.reduction = "integrated.rpca",
                          dims = 1:5, k.anchor = 20, reference = which(Layers(object, search = "data") %in% c("data.ZT0")),
                          verbose = F)
# cluster the integrated data
object <- FindNeighbors(object, reduction = "integrated.rpca", dims = 1:30)
object <- FindClusters(object, resolution = 0.1)
object <- RunUMAP(object, reduction = "integrated.rpca", dims = 1:30, return.model = T, verbose = F)
object[["sketch"]] <- JoinLayers(object[["sketch"]])


###PROJECT onto the rest of the cells
##this is now all the cells sequenced from his experiment! 
object[["sketch"]] <- split(object[["sketch"]], f = object$sample)

object <- ProjectIntegration(object = object, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")
object <- ProjectData(object = object, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.rpca.full",
                      full.reduction = "integrated.rpca.full", dims = 1:30)
object <- RunUMAP(object, reduction = "integrated.rpca.full", dims = 1:30, reduction.name = "umap.full",
                  reduction.key = "UMAP_full_")



#######LOUPER File for Manual Annotations#######

DefaultAssay(object) <- 'RNA'
object[['sketch']] <- NULL
object <- JoinLayers(object)

object <- DietSeurat(object, dimreducs = 'umap.full')
object 

create_loupe_from_seurat(object,
                         output_dir = getwd(),
                         output_name = 'IntegratedALLCELLS ZT0.6.12.18')



######Annotations File from Manual Annotations#####

DefaultAssay(object) <- 'RNA'
object[['sketch']] <- NULL
object <- JoinLayers(object)

annotations <- read.csv('integrated.rpca.full.res.1 v5.csv')
cells <- colnames(object)
table(cells == annotations$Barcode) #ALL TRUE 
table(annotations$integrated.rpca.full.res.1)
object$annotations <- annotations$integrated.rpca.full.res.1

table(object$annotations)
DimPlot(object)
