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

setwd("~/WGCNA_course/WGCNA/")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

nCPU = 8
enableWGCNAThreads(nThreads = nCPU)

### =============== Generate Edge- and Node-files for a set of genes of your choice, which can be used as Cytoscape input ================== ###
dir.create("./Data/cytoscape")

#load gene list with your genes of interest
genes <- read.table("./Data/gene_list.txt", header = F)

#load adjacency matrix of your network
#Takes a lot of RAM!
load(file = "Data/WGCNA_course_Net1_adjacency_matrix.RData")

#adjacency-treshold, genes that are connected with with edge-weights lower than this threshold will be excluded from the Node- and Edge-lists
adj_threshold = 0.05

#Extract genes of your choice from the adjacency matrix
genes <- genes[,1]
modadj = adjacency[genes, genes]
modadj = as.table(modadj)

#Export the network into edge and node list files that Cytoscape can read
cyt = exportNetworkToCytoscape(modadj,
                               edgeFile = "./Data/cytoscape/edges_subnetwork1.txt",
                               nodeFile = "./Data/cytoscape/nodes_subnetwork1.txt",
                               weighted = TRUE,
                               threshold = adj_threshold,
                               nodeNames = rownames(genes),
                               nodeAttr = NULL)

