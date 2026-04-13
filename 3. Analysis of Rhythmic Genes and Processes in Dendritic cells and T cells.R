###############################################################
###############################################################
#      R script corresponding to Maples et al.                #
# A single-cell atlas of intestinal immune cells across the   # 
#   day-night cycle reveals dynamic populations               #
###############################################################
###############################################################

##All packages used in this script can found cited in the Key Resource Table attached to the manuscript##

##########################THIRD SCRIPT OF FOUR#############################################
####### Analysis of circadian genes and processes in Dendritic cells and T cells ##########


########DC AND T CELL AGGREGATE EXPRESSION BY CELL TYPE##########

rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Matrix)
# set this option when analyzing large datasets
options(future.globals.maxSize = 300000 * 1024^2)
set.seed(1234)


#loading in annotated single cell object 
#see "1. Integration and Annotation of Circacdian Small Intestine Atlas.R" for parameters used
#and object details
set.seed(1234)
load('IntegratedALLCELLS.RData')

#######get aggregate expression file for Dendritic cells
dcs <- subset(object, subset = annotations %in% c('CD103+ cDC2s','pDCs', 'cDC1s'))

tmp <- AggregateExpression(dcs, group.by = 'time', return.seurat = TRUE)
write.table(tmp[['RNA']]$data, file = 'AggregateExpression Dendritic cells by Time.txt')

dcs$time_type <- paste(dcs$time, dcs$annotations, sep = '_')
tmp <- AggregateExpression(dcs, group.by = 'time_type', return.seurat = TRUE)
write.table(tmp[['RNA']]$data, file = 'AggregateExpression Dendritic cells by Time and Type.txt')

#######get aggregate expression file for T cells
t.list <- unique(object$annotations)[grep("T", unique(object$annotations))]
t.list <- t.list[-which(t.list == 'Transitional B-Cell')]
t.list <- c(t.list, 'NK-like')
t <- subset(object, subset = annotations %in% t.list)

t$time <- t$sample
table(t$time)

tmp <- AggregateExpression(t, group.by = 'time', return.seurat = TRUE)
write.table(tmp[['RNA']]$data, file = 'AggregateExpression T cells by Time.txt')
t$time_type <- paste(t$time, t$annotations, sep = '_')
tmp <- AggregateExpression(t, group.by = 'time_type', return.seurat = TRUE)
write.table(tmp[['RNA']]$data, file = 'AggregateExpression T cells by Time and Type.txt')


######Metacycle Analysis of Bulk T cell, Bulk DCs, and individual subtypes 

files <- list.files(path = "~/metainput", pattern = 'subtype.txt', all.files = FALSE, full.names = FALSE)
outD <- meta2d(infile=paste(files[i]), filestyle="txt", timepoints = c(0,6,12,18,24,30,36,42))


####################ssGSEA Analysis########################
#using same aggregate expression files generated above for ssGSEA pathway quantification
#using the R script Created: September 09, 2017 by Karsten Krug 
#please see ssgsea-gui.R for full functionality
#below are user defined parameters that were used. No modification was made to the base-code. 

## ssGSEA / PTM-SEA parameters
sample.norm.type    = "rank"              ## "rank", "log", "log.rank", "none" 
weight              = 0.75                ## value between 0 (no weighting) and 1 (actual data counts)
statistic           = "area.under.RES"    ## "Kolmogorov-Smirnov"
output.score.type   = "NES"               ## 'ES' or 'NES'
nperm               = 1e3                 ## No. of permutations
min.overlap         = 2                  ## minimal overlap between gene set and data
correl.type         = "z.score"           ## 'rank', 'z.score', 'symm.rank'
par                 = T                   ## use 'doParallel' package?
spare.cores         = 1                   ## No. of cores to leave idle
export.signat.gct   = T                   ## if TRUE gene set GCT files will be exported 
extended.output     = T                   ## if TRUE the GCT files will contain stats on gene set overlaps etc.   






#################CELLCHAT ANALYSIS###############

rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(CellChat)
options(future.globals.maxSize = 3000000 * 1024^2)
options(stringsAsFactors = FALSE)
set.seed(1234)


#loading in annotated single cell object 
#see "1. Integration and Annotation of Circacdian Small Intestine Atlas.R" for parameters used
#and object details
load('IntegratedALLCELLS.RData')


#isolate DC and T cells from the main Seurat object 
dc.t <- subset(object, subset = annotations %in% keep)
table(dc.t$annotations)
DefaultAssay(dc.t) <- "RNA"

dc.t <- SetIdent(dc.t, value = "annotations")
dc.t.0 <- subset(dc.t, subset = sample == 'ZT0')
dc.t.6 <- subset(dc.t, subset = sample == 'ZT6')
dc.t.12 <- subset(dc.t, subset = sample == 'ZT12')
dc.t.18 <- subset(dc.t, subset = sample == 'ZT18')

###now load all 4 time points as seperate cellchat objects
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
CellChatDB.use <- subsetDB(CellChatDB)

cellChat.zt0 <- createCellChat(object = dc.t.0, group.by = "annotations", assay = "RNA")
cellChat.zt6 <- createCellChat(object = dc.t.6, group.by = "annotations", assay = "RNA")
cellChat.zt12 <- createCellChat(object = dc.t.12, group.by = "annotations", assay = "RNA")
cellChat.zt18 <- createCellChat(object = dc.t.18, group.by = "annotations", assay = "RNA")

cellChat.zt0@DB <- CellChatDB.use
cellChat.zt0 <- subsetData(cellChat.zt0) # This step is necessary even if using the whole database
future::plan("multisession", workers = 10) # do parallel
cellChat.zt0 <- identifyOverExpressedGenes(cellChat.zt0)
cellChat.zt0 <- identifyOverExpressedInteractions(cellChat.zt0)
cellChat.zt0 <- computeCommunProb(cellChat.zt0, type = "triMean", population.size = TRUE) # this truncates at 25% by default
cellChat.zt0 <- filterCommunication(cellChat.zt0, min.cells = 10)
cellChat.zt0 <- computeCommunProbPathway(cellChat.zt0)
cellChat.zt0 <- aggregateNet(cellChat.zt0)

cellChat.zt6@DB <- CellChatDB.use
cellChat.zt6 <- subsetData(cellChat.zt6) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellChat.zt6 <- identifyOverExpressedGenes(cellChat.zt6)
cellChat.zt6 <- identifyOverExpressedInteractions(cellChat.zt6)
cellChat.zt6 <- computeCommunProb(cellChat.zt6, type = "triMean", population.size = TRUE) ##this truncates at 25% by default
cellChat.zt6 <- filterCommunication(cellChat.zt6, min.cells = 10)
cellChat.zt6 <- computeCommunProbPathway(cellChat.zt6)
cellChat.zt6 <- aggregateNet(cellChat.zt6)

cellChat.zt12@DB <- CellChatDB.use
cellChat.zt12 <- subsetData(cellChat.zt12) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellChat.zt12 <- identifyOverExpressedGenes(cellChat.zt12)
cellChat.zt12 <- identifyOverExpressedInteractions(cellChat.zt12)
cellChat.zt12 <- computeCommunProb(cellChat.zt12, type = "triMean", population.size = TRUE) ##this truncates at 25% by default
cellChat.zt12 <- filterCommunication(cellChat.zt12, min.cells = 10)
cellChat.zt12 <- computeCommunProbPathway(cellChat.zt12)
cellChat.zt12 <- aggregateNet(cellChat.zt12)

cellChat.zt18@DB <- CellChatDB.use
cellChat.zt18 <- subsetData(cellChat.zt18) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellChat.zt18 <- identifyOverExpressedGenes(cellChat.zt18)
cellChat.zt18 <- identifyOverExpressedInteractions(cellChat.zt18)
cellChat.zt18 <- computeCommunProb(cellChat.zt18, type = "triMean", population.size = TRUE) ##this truncates at 25% by default
cellChat.zt18 <- filterCommunication(cellChat.zt18, min.cells = 10)
cellChat.zt18 <- computeCommunProbPathway(cellChat.zt18)
cellChat.zt18 <- aggregateNet(cellChat.zt18)


cellChat.zt0 <- netAnalysis_computeCentrality(cellChat.zt0, slot.name = "netP") ## dataset 1
cellChat.zt6 <- netAnalysis_computeCentrality(cellChat.zt6, slot.name = "netP") ## dataset2
cellChat.zt12 <- netAnalysis_computeCentrality(cellChat.zt12, slot.name = "netP") ## dataset3
cellChat.zt18 <- netAnalysis_computeCentrality(cellChat.zt18, slot.name = "netP") ## dataset4
object.list <- list(ZT0 = cellChat.zt0, ZT6 = cellChat.zt6, ZT12 = cellChat.zt12, ZT18 = cellChat.zt18)
cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix=TRUE)








