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

### ----------------------------------------- STAGE IV -------------------------------------------------
### ---------------------------- Module- and gene-treatment correlations -------------------------------

#load expression data
load(file = "Data/dataInput.RData")

#Load network data
load(file = "Data/WGCNA_course_Net1_construction.RData")

#Load Module-Eigengenes (MEs) and order them according to similarities
MEs <- orderMEs(moduleEigengenes(data, moduleColors)$eigengenes)

## ------------------------------------------- STEP I --------------------------------------------------
##Generate an Eigengene-treatment table with the ME correlation values to each treatment
#MEs are artificial expression profiles representing the first principal component of each module

#From WGCNA_stage1.R "data" contains your samples (e.g. RNAseq libraries) as rows and genes as columns
#The order of your samples is preserved like in your initial input expression matrix
#In case you're unsure, check with rownames(data) manually

control = colMeans(MEs[1:2,])
cold = colMeans(MEs[3:4,])
heat = colMeans(MEs[5:6,])

MEs1 = data.frame(rbind(control, cold, heat))

write.table(MEs1, "Data/WGCNA_course_Net1_eigengene-treatment.txt", sep = "\t", quote = F)

## ------------------------------------------ STEP II --------------------------------------------------
##Generate a matrix with treatments for module-treatment correlation
#vector containing treatment names as many times as there are samples for each treatment
treatment <- c(rep("control", 2),
               rep("cold", 2),
               rep("heat", 2))

#dataframe with a column "treatment" containing the treatment names from before
treatment <- data.frame(treatment = as.factor(treatment))
#to this dataframe, add the sample names from the input data
rownames(treatment) <- rownames(data)

#create a modelmatrix, which states that control samples correlate maximally (1) with control treatment
#Mind the order of treatments in this matrix, for the renaming
treat <- model.matrix(~., data = treatment,
                      contrasts.arg = lapply(treatment, contrasts, contrasts = F))

#remove the first column of this modelmatrix, which is the "intercept"-column
treat <- treat[,-1]
#MEs matrix is correlated with treat modelmatrix
moduleTraitCor <- WGCNA::cor(MEs, treat, use = "p")
#Calculates student asymptotic p-value for correlations
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,length(treat))
#nicer column names
colnames(moduleTraitCor) <- c("cor_control", "cor_cold", "cor_heat")
#nicer row names
rownames(moduleTraitCor) <- rownames(moduleTraitCor) %>% substr(.,3,nchar(rownames(moduleTraitCor)))
#Do the same for the Pvalue table
colnames(moduleTraitPvalue) <- c("pval_control", "pval_cold", "pval_heat")
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
names(trait) = c("control", "cold", "heat")

##Write "gene to treatment correlation"- and "correlation pvalue"-files
dir.create("./Data/gene-treatment_cor")

col_names <- colnames(trait)

for (i in col_names) {
  
  #create a dataframe containing only one treatment column of choice
  selection = as.data.frame(trait[i])
  
  #nicer name for the chosen column
  names(selection) = i
  
  #correlate transcripts with the selected treatment
  geneTraitSignificance = as.data.frame(cor(data, selection, use = "p"))
  #get p-values for these correlations
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(data)))
  
  #nicer column titles
  names(geneTraitSignificance) = paste("cor.", names(selection), sep="")
  names(GSPvalue) = paste("pval.", names(selection), sep="")
  
  #save the produced tables for later merging in excel
  write.table(geneTraitSignificance, file = paste("./Data/gene-treatment_cor/", i, "_geneTraitSignificance.txt", sep = ""), sep = "\t", quote = F)
  write.table(GSPvalue, file = paste("./Data/gene-treatment_cor/", i, "_GSPvalue.txt", sep = ""), sep = "\t", quote = F)
  
}

#merge the individual files for convenience by
#listing all files
mult_files = list.files("./Data/gene-treatment_cor/", pattern="*[Significance,GSPvalue].txt")
mult_files <- paste("./Data/gene-treatment_cor/", mult_files[1:length(mult_files)], sep = "")
#reading all files
myfilelist <- lapply(mult_files, read.table)
#assigning variable names to the files
names(myfilelist) <- list.files("./Data/gene-treatment_cor/", pattern="*[Significance,GSPvalue].txt", full.names=FALSE)
#rbind all the files together
gene_treatment_cor <- as.data.frame(lapply(myfilelist, rbind))
#save the file
write.table(gene_treatment_cor, file = "./Data/Combined_gene-treatment_cor.txt", sep = "\t", quote = F)
#after combining the files, remove the directory with the individual files
unlink("./Data/gene-treatment_cor", recursive = TRUE)

## ------------------------------------------- STEP IV ------------------------------------------------
#Extract top10 correlated genes for each module, based on the correlation score of a gene to a module-eigengene

#vector of nice module names
modNames = substring(names(MEs), 3)
#membership of a gene to a module, expressed in pearson-correlation
geneModuleMembership = as.data.frame(cor(data, MEs, use = "p"))
#pvalues of these correlations
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(data)))
#changing column titles from ME to MM
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
#combine correlation and pvalue tables
geneModuleMembership <- cbind(geneModuleMembership, MMPvalue)
#save the produced tables
write.table(geneModuleMembership, "./Data/geneModuleMemebership.txt", sep = "\t", quote = F)

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

write.table(hubs, "./Data/WGCNA_course_Net1_Top_10_ModuleMembership_genes.txt", sep = "\t", quote = F)

## ------------------------------------------- STEP V ------------------------------------------------
#Extract the top1 HUB-gene for each module, based on the connectivity score of a gene

hub = chooseTopHubInEachModule(data, moduleColors, omitColors = "grey",
                               power = 1, type = "unsigned")

write.table(hub,"./Data/Module_top1_HUB-genes.txt", quote = F, col.names = F)
