###############################################################
###############################################################
#      R script corresponding to Maples et al.                #
# A single-cell atlas of intestinal immune cells across the   # 
#   day-night cycle reveals dynamic populations               #
###############################################################
###############################################################

##All packages used in this script can found cited in the Key Resource Table attached to the manuscript##

################################SECOND SCRIPT OF FOUR###########################################
#######Overlap of Xu et al Reference Small Intestine Dataset with newly generated data##########



rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Matrix)
# set this option when analyzing large datasets
options(future.globals.maxSize = 300000 * 1024^2)
set.seed(1234)

########METADATA XU ET AL (object name: clusters)#################
##load in the cell count matrices and meta data from Xu et al GEO submission, GSE124880
barcodes <- read.csv(gzfile('GSE124880_PP_LP_mm10_count_barcodes.tsv.gz'), header = FALSE)
clusters <- read.csv('Barcode_cluster.csv')

clusters$CellID_shortened <- ''
##to match them we need to take off the first 16 characters from the clusters file 
for (i in 1:length(clusters$CellID)){
  clusters$CellID_shortened[i] <- substr(clusters$CellID[i], 17, nchar(clusters$CellID[i]))
}

table(barcodes$V1 == clusters$CellID_shortened) #awesome! they all match now 


#location ids 
lp <- read.table(gzfile('/project/immunology/Hooper_lab/s189983/scRNAseq_Bobby/Integration with xu et al/Food-Allergy-LP_normalized_matrix.txt.gz'), header = TRUE)
pp <- read.table(gzfile('/project/immunology/Hooper_lab/s189983/scRNAseq_Bobby/Integration with xu et al/Food-Allergy-PP_normalized_matrix.txt.gz'), header = TRUE)

##add these to the clusters dataframe 
clusters$dottedname <- ''
for (i in 1:length(clusters$CellID)){
  x <- strsplit(clusters$CellID[i],'-', fixed = TRUE)[[1]]
  dottedname <- x[1]
  for (l in 2:length(x)){
    dottedname <- paste0(dottedname, ".", x[l])
  }
  clusters$dottedname[i] <- dottedname
}

table(clusters$dottedname %in% colnames(pp))
table(clusters$dottedname %in% colnames(lp))




clusters$location <- ''
for (i in 1:length(clusters$dottedname)){
  if (clusters$dottedname[i] %in% colnames(lp)){
    clusters$location[i] <- 'LP'
  }
  if (clusters$dottedname[i] %in% colnames(pp)){
    clusters$location[i] <- 'PP'
  }
}
table(clusters$location)


#rename the clusters based on their cell identity from Supplemental Figure they have in their paper. 
clusternumber <- clusters$cluster
clusternumber[which(clusternumber == '1')] <- 'Resting CD4+ T cell'
clusternumber[which(clusternumber == '5')] <- 'CD8+ T cell'
clusternumber[which(clusternumber == '7')] <- 'Activated CD4+ T cell'
clusternumber[which(clusternumber == '13')] <- 'GREY CD4+ T cell (low UMI count)'
clusternumber[which(clusternumber == '17')] <- 'gd T cell (Xcl1+)'
clusternumber[which(clusternumber == '18')] <- 'GREY Unresolved'
clusternumber[which(clusternumber == '19')] <- 'NKT cell'
clusternumber[which(clusternumber == '38')] <- 'gd T cell (Gzma+)'
clusternumber[which(clusternumber == '42')] <- 'GREY Doublets'
clusternumber[which(clusternumber == '45')] <- 'GREY T precursor-like cell'
clusternumber[which(clusternumber == '2')] <- 'Resting B cell'
clusternumber[which(clusternumber == '9')] <- 'GC B cell (DZ)'
clusternumber[which(clusternumber == '10')] <- 'GC B cell (LZ)'
clusternumber[which(clusternumber == '20')] <- 'Plasma cell'
clusternumber[which(clusternumber == '23')] <- 'GREY Resting B cell (low UMI count)'
clusternumber[which(clusternumber == '29')] <- 'GREY Doublets'
clusternumber[which(clusternumber == '43')] <- 'GREY Doublets'
clusternumber[which(clusternumber == '3')] <- 'ILC3'
clusternumber[which(clusternumber == '4')] <- 'LTi cell'
clusternumber[which(clusternumber == '8')] <- 'NK cell'
clusternumber[which(clusternumber == '11')] <- 'ILC1'
clusternumber[which(clusternumber == '12')] <- 'ILC2'
clusternumber[which(clusternumber == '15')] <- 'GREY ILC3 (low UMI count)'
clusternumber[which(clusternumber == '16')] <- 'GREY LTi cell (low UMI count)'
clusternumber[which(clusternumber == '27')] <- 'GREY NK cell (low UMI count)'
clusternumber[which(clusternumber == '37')] <- 'GREY Unresolved'
clusternumber[which(clusternumber == '6')] <- 'pDC'
clusternumber[which(clusternumber == '14')] <- 'DC (CD103+CD11b+)'
clusternumber[which(clusternumber == '21')] <- 'DC (CD103+CD11b-)'
clusternumber[which(clusternumber == '24')] <- 'DC (CD103- C1)'
clusternumber[which(clusternumber == '31')] <- 'DC (CD103- C2)'
clusternumber[which(clusternumber == '34')] <- 'GREY Unresolved'
clusternumber[which(clusternumber == '41')] <- 'GREY Doublets'
clusternumber[which(clusternumber == '22')] <- 'Macrophage'
clusternumber[which(clusternumber == '26')] <- 'Mast cell'
clusternumber[which(clusternumber == '28')] <- 'GREY GC B cell (low UMI count)'
clusternumber[which(clusternumber == '36')] <- 'Neutrophil'
clusternumber[which(clusternumber == '40')] <- 'GREY Macrophage (low UMI count)'
clusternumber[which(clusternumber == '44')] <- 'Basophil'
clusternumber[which(clusternumber == '25')] <- 'Endothelial cell'
clusternumber[which(clusternumber == '32')] <- 'Stromal cell (DN)'
clusternumber[which(clusternumber == '33')] <- 'Fibroblast'
clusternumber[which(clusternumber == '46')] <- 'GREY Lymphatic endothelial-like cell'
clusternumber[which(clusternumber == '30')] <- 'GREY Epithelial cell C1'
clusternumber[which(clusternumber == '35')] <- 'GREY Epithelial cell C2'
clusternumber[which(clusternumber == '39')] <- 'GREY Epithelial cell (low UMI count)'
clusters$cellidentity <- clusternumber







########################COUNT MATRICES XU ET AL LOADED IN########################

#downloaded count matrices from GEO submission GSE124880
barcodes <- read.csv(gzfile("GSE124880_PP_LP_mm10_count_barcodes.tsv.gz"), header = FALSE)
features <- read.csv(gzfile("GSE124880_PP_LP_mm10_count_genes.tsv.gz"), header = FALSE)
matrix <- Matrix::readMM("GSE124880_PP_LP_mm10_count_matrix.mtx.gz")
colnames(matrix) <- barcodes$V1
rownames(matrix) <- features$V1

xu <- CreateSeuratObject(counts = matrix, min.cells = 0, min.features = 0) #since uploaded GSE data is already filtered according to authors

xu <- NormalizeData(xu)
xu <- FindVariableFeatures(xu)
xu <- ScaleData(xu, features = rownames(xu))
xu <- RunPCA(xu)
xu <- FindNeighbors(xu, dims = 1:30)
xu <- FindClusters(xu)
xu <- RunUMAP(xu, dims = 1:30, return.model = TRUE)

table(colnames(xu) == clusters$CellID_shortened)
xu$cellidentity <- clusters$cellidentity
xu$location <- clusters$location
xu$published_clusters <- clusters$cluster



###############LOAD IN MANUSCRIPTS DATASET (object name: object, created in script #1)###############

#loading in annotated single cell object 
#see "1. Integration and Annotation of Circacdian Small Intestine Atlas.R" for parameters used. 

set.seed(1234)
load('IntegratedALLCELLS.RData')

query <- DietSeurat(object, assays = 'RNA',layers = c('counts','data'), dimreducs = c('integrated.rpca.full','umap.full'))
# in order to integrate dimension reduction plots need to the same naming convention 
names(query@reductions)[names(query@reductions) == "integrated.rpca.full"] <- "pca"
names(query@reductions)[names(query@reductions) == "umap.full"] <- "umap"

#######combine both datasets and integrate on the same UMAP 
combined <- merge(xu, query, add.cell.ids = c('Ref','Query'))
barcodes <- colnames(combined)
barcodes[grep('Ref', barcodes)] <- 'Reference'
barcodes[grep('Query', barcodes)] <- 'Query'
table(barcodes)
combined$projects <- barcodes



# run standard analysis workflow
combined <- NormalizeData(combined)
combined[['RNA']] <- JoinLayers(combined[['RNA']])
combined[['RNA']] <- split(combined[['RNA']], f=combined$projects)
combined <- FindVariableFeatures(combined, verbose = FALSE)

combined <- SketchData(object= combined, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")

DefaultAssay(combined) <- "sketch"
combined <- FindVariableFeatures(combined, verbose = F)
combined <- ScaleData(combined, verbose = F)
combined <- RunPCA(combined, verbose = F)
# integrate the datasets
combined <- IntegrateLayers(combined, method = RPCAIntegration, orig = "pca", new.reduction = "integrated.rpca",
                            dims = 1:2, k.anchor = 20, reference = which(Layers(combined, search = "data") %in% c("data.Reference")),
                            verbose = F)

# cluster the integrated data
combined <- FindNeighbors(combined, reduction = "integrated.rpca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.1)
combined <- RunUMAP(combined, reduction = "integrated.rpca", dims = 1:30, return.model = T, verbose = F)
combined[["sketch"]] <- JoinLayers(combined[["sketch"]])


###PROJECT onto the rest of the cells
##this is now all the cells sequenced from this manuscript
combined[["sketch"]] <- split(combined[["sketch"]], f = combined$projects)

combined <- ProjectIntegration(object = combined, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")
combined <- ProjectData(object = combined, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.rpca.full",
                        full.reduction = "integrated.rpca.full", dims = 1:30)
combined <- RunUMAP(combined, reduction = "integrated.rpca.full", dims = 1:30, reduction.name = "umap.full",
                    reduction.key = "UMAP_full_")
