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

setwd("/media/Data2/marc2/Biscutella_genome-based_RNAseq/WGCNA")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

nCPU = 12
enableWGCNAThreads(nThreads = nCPU)


### ----------------------------------------- STAGE IV -------------------------------------------------
### ------------------------------------------ Tables --------------------------------------------------

#load expression data
load(file = "Data/dataInput.RData")

# Load network data
load(file = "Data/2ndSE_GBGE_Net2_construction.RData")

#Define numbers of genes and samples
nGenes = ncol(data)
nSamples = nrow(data)

#Load eigengenes
{
  MEs <- moduleEigengenes(data, moduleColors)$eigengenes #calculates module eigengenes
  MEs0 <- moduleEigengenes(data, moduleColors)$eigengenes
  MEs <- orderMEs(MEs0)
}

## ------------------------------------------- STEP I --------------------------------------------------
##Generate an Eigengene-treatment table, Eigengenes are defined as the first principal component of each module
#Average the correlation of the eigengenes for each module across all replicates of a treatment

#make sure to have the correct row-order of treatments here!
{
  cold = colMeans(MEs[1:4,])
  control = colMeans(MEs[5:8,])
  drought = colMeans(MEs[9:11,])
  heat = colMeans(MEs[12:14,])
  herbivory = colMeans(MEs[15:17,])
}

MEs1 = data.frame(rbind(cold, control, drought, heat, herbivory))

write.table(MEs1, "Data/2ndSE_GBGE_Net2_eigengene-treatment.txt", sep = "\t", quote = F)

#remove the vectors containing eigengene averages from the environment
rm(cold, control, drought, heat, herbivory)

## ------------------------------------------ STEP II --------------------------------------------------
##Generate a matrix with treatments for module-treatment correlation
#vector containing treatment names as many times as there are samples for each treatment
treatment <- c(rep("cold", 4),
               rep("control", 4),
               rep("drought", 3),
               rep("heat", 3),
               rep("herbivory", 3))

#dataframe with a column "treatment" containing the treatment names from before
treatment <- data.frame(treatment = as.factor(treatment))
#to this dataframe, add the sample names from the input data
rownames(treatment) <- rownames(data)

#create a modelmatrix, which states that control samples correlate maximally (1) with control treatment
#Mind the order of treatments in this matrix, for the renaming
treat <- model.matrix(~., data = treatment,
                      contrasts.arg = lapply(treatment, contrasts, contrasts = F))

#remove the first column of this modelmatrix, which is the apparently useless "intercept"-column
treat <- treat[,-1]
#MEs matrix is correlated with treat modelmatrix
moduleTraitCor <- WGCNA::cor(MEs, treat, use = "p")
#Calculates student asymptotic p-value for correlations
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,length(treat))
#nicer column names
colnames(moduleTraitCor) <- c("cor_cold", "cor_control", "cor_drought", "cor_heat", "cor_herbivory")
#nicer row names
rownames(moduleTraitCor) <- rownames(moduleTraitCor) %>% substr(.,3,nchar(rownames(moduleTraitCor)))
#Do the same for the Pvalue table
colnames(moduleTraitPvalue) <- c("pval_cold", "pval_control", "pval_drought", "pval_heat", "pval_herbivory")
rownames(moduleTraitPvalue) <- rownames(moduleTraitPvalue) %>% substr(.,3,nchar(rownames(moduleTraitPvalue)))
#Combine correlation and pvalue tables
module_treatment_cor <- cbind(moduleTraitCor, moduleTraitPvalue)
#save data for STAGE V and write a table
save(MEs, moduleTraitCor, moduleTraitPvalue, file = "Data/module_trait_cor.RData")
write.table(module_treatment_cor, "Data/Module-treatment_cor.txt", sep = "\t", quote = F)

## ------------------------------------------ STEP III -------------------------------------------------
##Correlate each transcript to each treatment
#create a treatment dataframe from the modelmatrix
trait=as.data.frame(treat)
#nicer column names, ordered as in "treat"
names(trait) = c("cold", "control", "drought", "heat", "herbivory")

##Write "gene to treatment correlation"- and "correlation pvalue"-files
dir.create("gene-treatment_cor")

col_names <- colnames(trait)

for (i in col_names) {

  #create a dataframe containing only one treatment column of choice
  selection = as.data.frame(trait[i])

  #nicer name for the chosen column
  names(selection) = i

  #correlate transcripts with the selected treatment
  geneTraitSignificance = as.data.frame(cor(data, selection, use = "p"))
  #get p-values for these correlations
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

  #nicer column titles
  names(geneTraitSignificance) = paste("cor.", names(selection), sep="")
  names(GSPvalue) = paste("pval.", names(selection), sep="")

  #save the produced tables for later merging in excel
  write.table(geneTraitSignificance, file = paste("gene-treatment_cor/", i, "_geneTraitSignificance.txt", sep = ""), sep = "\t", quote = F)
  write.table(GSPvalue, file = paste("gene-treatment_cor/", i, "_GSPvalue.txt", sep = ""), sep = "\t", quote = F)

}

#merge the individual files for convenience by
#listing all files
mult_files = list.files("gene-treatment_cor/", pattern="*[Significance,GSPvalue].txt")
mult_files <- paste("gene-treatment_cor/", mult_files[1:length(mult_files)], sep = "")
#reading all files
myfilelist <- lapply(mult_files, read.table)
#assigning variable names to the files
names(myfilelist) <- list.files("gene-treatment_cor/", pattern="*[Significance,GSPvalue].txt", full.names=FALSE)
#rbind all the files together
gene_treatment_cor <- as.data.frame(lapply(myfilelist, rbind))
#save the file
write.table(gene_treatment_cor, file = "Data/Combined_gene-treatment_cor.txt", sep = "\t", quote = F)
#after combining the files, remove the directory with the individual files
unlink("gene-treatment_cor", recursive = TRUE)


## ------------------------------------------- STEP IV ------------------------------------------------
#Extract top10 correlated genes for each module, based on the correlation score of a gene to a module-eigengene

#vector of nice module names
modNames = substring(names(MEs), 3)
#membership of a gene to a module, expressed in correlation
geneModuleMembership = as.data.frame(cor(data, MEs0, use = "p"))
#pvalues of these correlations
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
#changing column titles from ME to MM
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
#combine correlation and pvalue tables
geneModuleMembership <- cbind(geneModuleMembership, MMPvalue)
#save the produced tables
write.table(geneModuleMembership, "Data/geneModuleMemebership.txt", sep = "\t", quote = F)

#Extract HUB-genes
clr = sort(unique(moduleColors))
hubs = character()

for (i in 1:length(clr)) {
  dt = geneModuleMembership[,i]
  kk = order(dt, decreasing = T)
  nm = rownames(geneModuleMembership)[kk]
  top10 = nm[1:10] #chose a number of top HUB-genes to extract from each module
  hubs = cbind(hubs,top10)
}

hubs = as.data.frame(hubs)
colnames(hubs) = clr

write.table(hubs,"Data/2ndSE_GBGE_Net2_Top_10_ModuleMembership_genes.txt", sep = "\t", quote = F)

## ------------------------------------------- STEP V ------------------------------------------------
#Extract the top1 HUB-gene for each module, based on the connectivity score of a gene

hub = chooseTopHubInEachModule(data, moduleColors, omitColors = "grey",
                               power = 10, type = "unsigned")

write.table(hub,"Data/Module_top1_HUB-genes.txt", quote = F, col.names = F)
