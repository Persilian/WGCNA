### ----------------------------- Setup Workspace ------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", quietly = TRUE)

#BiocManager::install("tidyverse")
#BiocManager::install("WGCNA")

library(WGCNA)
library(tidyverse)


setwd("/media/Data2/marc2/Biscutella_genome-based_RNAseq/WGCNA/Net2/DEG_subnetwork")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

nCPU = 12
enableWGCNAThreads(nThreads = nCPU)

### ============================ create edge and node files for Cytoscape ============================ ###
#requires an adjacency matrix (can be extracted from the main network using exa_v2.MATRIX.sh)
#optionally one can provide a data.frame with node-attributes

#load adjacency matrix of the subnetwork
adj <- read.table("blaev_all_DEG_transcripts_adjacency.txt", header = TRUE)

probes <- colnames(adj)

#load node attribute data
node_att <- read.table("blaev_DEG_subnetwork_node_attributes.txt", header = TRUE)


exportNetworkToCytoscape(
   adj,
   edgeFile = "cytoscape/blaev_GBGE_all_DEGs_thr01_edges.txt", #thr05 means threshold = 0.5
   nodeFile = "cytoscape/blaev_GBGE_all_DEGs_thr01_nodes.txt",
   weighted = TRUE,
   threshold = 0.1, #how to set the threshold? Let's leave it at 0.5 for now, if too many connections, increase it and vice-versa
   nodeNames = probes,
   nodeAttr = node_att)


