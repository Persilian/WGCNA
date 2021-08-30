### ----------------------------- Setup Workspace ------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", quietly = TRUE)

#BiocManager::install("tidyverse")
#BiocManager::install("WGCNA")
#BiocManager::install("ggdendro")
#BiocManager::install("gplots")
#BiocManager::install("reshape2")
#BiocManager::install("patchwork")
#BiocManager::install("circlize")

library(WGCNA)
library(tidyverse)
library(ggdendro)
library(gplots)
library(reshape2)
library(patchwork)
library(circlize)


setwd("~/WGCNA")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

nCPU = 8
enableWGCNAThreads(nThreads = nCPU)

### ----------------------------------------- STAGE V -------------------------------------------------
### ----------------------------------------- Heatmaps ------------------------------------------------


## ------------------------------------------- STEP I -------------------------------------------------
##Module-treatment correlation heatmap
#load data
load(file = "Data/module_trait_cor.RData")
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

#define color gradient
col_fun = colorRamp2(c(-100, 0, 100), c("#648FFF", "white", "#FFB000"))
heatmappalette <- col_fun(seq(-150, 150))

#Open a large plot window
sizeGrWindow(12,20)

pdf("Plots/detailed_module-treatment-heatmap.pdf") #open pdf to contain the plot
par(mar = c(6, 8.5, 3, 3)) #define graph size parameters

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = substr(colnames(moduleTraitCor), 5, nchar(colnames(moduleTraitCor))),
               yLabels = names(MEs),
               ySymbols = substring(names(MEs), 3),
               colorLabels = FALSE,
               colors = heatmappalette,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.2,
               cex.lab.y = 0.5,
               zlim = c(-1,1),
               main = paste("Module-Treatment correlations"))

dev.off() #finish and save the plot opened above

## ------------------------------------------- STEP II ------------------------------------------------
##SIMPLE_Module-treatment correlation heatmap for easier choice of interesting modules
#Modules with a correlation of higher than "cor_threshold" and a p-value of lower than "p_threshold" are colored differently

#Matrix of logicals used for the simplified heatmap
cor_treshold <- 0.7 #chose correlation threshold to be depicted
p_threshold <- 0.05 #chose p-avlue threshold to be depicted
cor <- abs(moduleTraitCor) >= cor_treshold
p <- moduleTraitPvalue < p_threshold
xx <- cor & p
mt <- as.matrix(xx)
dt <- as.integer(mt)
dt <- matrix(dt, nrow = dim(moduleTraitCor)[1], ncol = dim(moduleTraitCor)[2])
colnames(dt) <- colnames(mt)
rownames(dt) <- rownames(mt)

#Open a large plot window
sizeGrWindow(12,20)

pdf("Plots/simple_module-treatment-heatmap.pdf") #open pdf to contain the plot
par(mar = c(6, 8.5, 3, 3)) #define graph size parameters

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = dt,
               xLabels = substr(colnames(dt), 5, nchar(colnames(dt))),
               yLabels = names(MEs),
               ySymbols = substring(names(MEs), 3),
               colorLabels = FALSE,
               colors = heatmappalette,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.2,
               cex.lab.y = 0.5,
               zlim = c(-1,1),
               main = paste("Module-Treatment correlations of R^2 >= ", cor_treshold))
dev.off()
