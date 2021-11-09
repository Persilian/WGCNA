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

setwd("~/WGCNA")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

nCPU = 8
enableWGCNAThreads(nThreads = nCPU)

### ------------------------------------------ STAGE III -----------------------------------------------
### ---------------------------- Adjacency matrix and gene connectivity --------------------------------

## -------------------------------------------- STEP I -------------------------------------------------
#Generate the adjacency matrix of your previously built network, for later extractions of sub-networks.

#load expression data from STAGE I
load(file = "Data/dataInput.RData")

##calculate adjacency matrix
#Based on the soft-thresholding choice in STAGE I, we choose the same power again.
#CAREFUL - DEFAULT of adjacency() is "unsigned"

adjacency = adjacency(data, type = "unsigned", power = 1)
save(adjacency, file = "Data/WGCNA_course_Net1_adjacency_matrix.RData")

## ---------------------------------------------- STEP II ----------------------------------------------
#Compute gene connectivity for each gene/isoform in the network. Requires adjacency matrix. 

#load adjacency matrix if you run this step without creating the adjacency matrix
#load(file = "Data/WGCNA_course_Net1_adjacency_matrix.RData")

#load network construction data created in WGCNA_stage2.R
load(file = "Data/WGCNA_course_Net1_construction.RData")

#The function intramodularConnectivity computes the whole network connectivity kTotal,
#the within module connectivity kWithin, kOut=kTotal-kWithin, and kDiff=kIn-kOut=2*kIN-kTotal
connectivity <- intramodularConnectivity(adjacency, net$colors)

write.table(connectivity, "Data/WGCNA_course_Net1_gene_connectivity.txt", sep = "\t", quote = F)
