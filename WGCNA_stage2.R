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

### ------------------------------------------ STAGE II -----------------------------------------------
### ------------------------------------ Network construction -----------------------------------------

##--------------------------------------------- STEP I ------------------------------------------------
##---------------------------------- One-step Network construction ------------------------------------

#load expression data from STAGE I
load(file = "Data/dataInput.RData")

#Based on the soft-thresholding choice in STAGE I, we chose an adequate power.
#power = 1

#We use "pearson" correlation (linear correlations) as this is a more intuitive and conservative approach
#networkType = "signed" means that strongly negative correlated gene expression profiles will not be connected
#For RNAseq data it makes biological sense to generate "unsigned" networks, where negative correlated expression profiles will be connected
#We would like to have as few "Blocks" as possible, generally you will have around 50'000 genes in your expression data from stage 1.
#Depending on how much RAM you have available, set maxBlockSize to something around 10'000.

#One-step network construction and module detection
net = blockwiseModules(data, power = 1, corType = "pearson", networkType = "unsigned",
                        nThreads = nCPU,
                        maxBlockSize = 10000,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.1,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        #loadTOM = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "Data/TOM.mergecutheight.0.1.WGCNA_course_Net1",
                        verbose = 3)

#number of modules and module sizes
sum = table(net$colors)
write.table(sum,"Data/Gene_count.txt", quote = F, row.names = F)

#saving the environment and remove trash module "grey"
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs0 = moduleEigengenes(data, moduleColors)$eigengenes
MEs = removeGreyME(MEs0,  greyMEName = paste(moduleColor.getMEprefix(), "grey", sep=""))
geneTree = net$dendrograms[[1]]
save(net, MEs, moduleLabels, moduleColors, geneTree,
     file = "Data/WGCNA_course_Net1_construction.RData")

##--------------------------------------------- STEP II -----------------------------------------------
##-------------------------------------- Network vizualization ----------------------------------------

#load data in case you have already computed the network
#load(file = "Data/WGCNA_course_Net1_construction.RData")

#plot the module dendrogram (network of genes)
#Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

#Plot the dendrogram and the module colors underneath
pdf(file = "Plots/Module_dendrogram.pdf", width = 12, height = 9)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                     "Module colors",
                     dendroLabels = FALSE, hang = 0.03,
                     addGuide = TRUE, guideHang = 0.05)
dev.off()

#plot the Module-Eigengene dendrogram (network of eigengenes)
#and the relationships among the eigengenes and the trait
pdf(file = "Plots/Module-Eigengene_dendrogram.pdf", width = 12, height = 9)
par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2),
                      cex.lab = 0.8,
                      xLabelsAngle = 90)
dev.off()

##------------------------------------------- STEP III ----------------------------------------------
##------------------------------ extract gene-lists from each module --------------------------------
#gene lists have not undergone cytoscape-extraction-thresholding at this point (full gene lists)
dir.create("Data/Gene_lists")

mod.clr = unique(moduleColors)
mod.clr = sort(mod.clr)
probes = colnames(data)

#Select modules
for (modules in mod.clr) {
  #Select module probes
  inModule = is.finite(match(moduleColors, modules))
  modProbes = probes[inModule]
  fileName = paste("Data/Gene_lists/gene_list_", modules, ".txt", sep="")
  write.table(as.data.frame(modProbes), file = fileName,
              row.names = FALSE, col.names = FALSE, quote = F)
}
