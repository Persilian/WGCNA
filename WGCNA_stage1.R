### ----------------------------- Load packages ------------------------------------

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

#Set your working directory
setwd("~/WGCNA")

#The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

nCPU = 12
enableWGCNAThreads(nThreads = nCPU)

#Create directories
dir.create("./Data")
dir.create("./Plots")

#Load expression data, for example a TMM-normalized TPM expression matrix from DESeq2 output
#Note that you should rename the .matrix into a .txt file and add a column-header "transcripts" to the first column
data <- read.table("WGCNA_course.isoform.TMM.EXPR.txt", header = TRUE)
data <- t(data)

### ----------------------------- STAGE I --------------------------------------------
### ----- Data Input and Preprocessing; Sample clustering; Soft-threshold choice -----

##------------------------------ STEP I ----------------------------------------------
##----------------------- Data Input and Preprocessing -------------------------------

##QC - checking for missing values
gsg = goodSamplesGenes(data, verbose = 2)
gsg$allOK

#If there are erroneous genes, remove them
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(data)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(data)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  data = data[gsg$goodSamples, gsg$goodGenes]
}

collectGarbage()

#saving data
save(data, file = "Data/dataInput.RData")

##--------------------------------- STEP II ---------------------------------------
##-------------------- Sample clustering and Soft-threshold choice ------------------

##Cluster samples
sampleTree = hclust(dist(data), method = "average")
# Plot the sample tree as a dendrogram
library(grDevices)
pdf(file = "Plots/1-sampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

##Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=100, by=2))

##Call the network topology analysis function (Takes time, so run all the commented parts only the first time)
sft = pickSoftThreshold(data, dataIsExpr = TRUE, powerVector = powers, verbose = 5)
write.table(sft, file = "Data/sft.txt", sep = "\t", quote = FALSE)
save(sft, file = "Data/sft.RData")
#load(file = "Data/sft.RData")

##Plot the results:
pdf(file = "Plots/2-thresholding.pdf", width = 12, height = 9)
par(mfrow = c(1,2))
cex1 = 0.9

##Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
      type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
      labels=powers,cex=cex1,col="red")

##this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

##Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
      main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

