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
### ---------------------------- Adjacency matrix generation (optional) --------------------------------

## -------------------------------------------- STEP I -------------------------------------------------
#For subsequent use with exa.MATRIX.sh and EGAD network analysis, the network adjacency matrix is produced here.

#load expression data from STAGE I
load(file = "Data/dataInput.RData")

##calculate adjacency matrix
#Based on the soft-thresholding choice in STAGE I, we chose the same power again.
#CAREFUL - DEFAULT of adjacency() is "unsigned"

adjacency = adjacency(data, type = "unsigned", power = 1)
save(adjacency, file = "Data/WGCNA_course_Net1_adjacency_matrix.RData")
#load(file = "Data/WGCNA_course_Net1_matrix.RData")
write.table(adjacency, "Data/WGCNA_course_Net1_adjacency_matrix.txt", sep = "\t", quote = F)


## ---------------------------------------------- STEP II ----------------------------------------------
#Compute gene connectivity. For each gene/isoform in the network, the total connections, 
#within module and outside module connections are calculated. Requires adjacency matrix. 

#load adjacency matrix and network
#load(file = "Data/WGCNA_course_Net1_adjacency_matrix.RData")
load(file = "Data/WGCNA_course_Net1_construction.RData")

#The function intramodularConnectivity computes the whole network connectivitykTotal,
#the within module connectivity kWithin, kOut=kTotal-kWithin, andkDiff=kIn-kOut=2*kIN-kTotal
connectivity <- intramodularConnectivity(adjacency, net$colors)

write.table(connectivity, "Data/WGCNA_course_Net1_gene_connectivity.txt", sep = "\t", quote = F)
