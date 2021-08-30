### ----------------------------- Setup Workspace ------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", quietly = TRUE)

#BiocManager::install("tidyverse")
#BiocManager::install("WGCNA")
#BiocManager::install("ggdendro")
#BiocManager::install("gplots")
#BiocManager::install("reshape2")
#BiocManager::install("patchwork")

library(WGCNA)
library(tidyverse)
library(ggdendro)
library(gplots)
library(reshape2)
library(patchwork)

setwd("~/WGCNA/")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

nCPU = 8
enableWGCNAThreads(nThreads = nCPU)

### ---------------------------------------- STAGE VI -------------------------------------------------
### ----------------------------- Cytoscape Node- and Edge-lists --------------------------------------

#load data
load(file = "Data/dataInput.RData")
load(file = "Data/WGCNA_course_Net1_construction.RData")

###parameters
#adjacency-treshold, genes that are connected with with edge-weights lower than this threshold will be excluded from the Node- and Edge-lists
adj_threshold = 0.05

## ------------------------------------------- STEP I -------------------------------------------------
#Generate Edge- and Node-files for each module, which can be used as Cytoscape input

load(file = "Data/WGCNA_course_Net1_adjacency_matrix.RData") #TAKES A LOT OF RAM

dir.create("Data/cytoscape")

clr = unique(moduleColors)
clr = sort(clr)
probes = colnames(data)

#Select modules
for (modules in clr) {
  
  #Select module probes
  inModule = is.finite(match(moduleColors, modules))
  modProbes = probes[inModule]
  
  #Select the corresponding Topological Overlap
  modadj = adjacency[inModule, inModule]
  dimnames(modadj) = list(modProbes, modProbes)
  modadj = as.table(modadj)

  #extract the adjacency submatrices of each module
  #write.table(modadj, file = paste("adjacency_", modules, ".txt", sep = ""))
  
  #Export the network into edge and node list files that Cytoscape can read
  cyt = exportNetworkToCytoscape(modadj,
                                 edgeFile = paste("Data/cytoscape/edges_", modules, ".txt", sep=""),
                                 nodeFile = paste("Data/cytoscape/nodes_", modules, ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = adj_threshold,
                                 nodeNames = modProbes,
                                 nodeAttr = moduleColors[inModule])
}

## ------------------------------------------- STEP II --------------------------------------------------
#Generate gene lists of modules after adjacency-thresholding

dir.create("Data/Thresholded_gene_lists")

for (i in 1: length(clr)) {
  modules = clr[i]
  path = paste("Data/cytoscape/edges_", modules, ".txt", sep = "")
  x = read.delim(path)
  dt = unique(c(as.character(x$fromNode,x$toNode)))
  fileName = paste("Data/Thresholded_gene_lists/gene_list_", modules, ".txt", sep = "")
  write.table(as.data.frame(dt), file = fileName,
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
