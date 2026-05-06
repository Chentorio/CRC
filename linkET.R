# Define input file for CIBERSORT
inputFile = "normalized_data_matrix.txt"      
# setwd("./") 

# Load CIBERSORT script (Ensure CIBERSORT.R is in the working directory)
source("CIBERSORT.R")       

# Execute CIBERSORT to estimate immune cell fractions
# 'ref.txt' is the signature matrix (e.g., LM22)
outTab = CIBERSORT("ref.txt", inputFile, perm=1000)

# Filter samples with P-value < 0.05 and select immune cell fraction columns
outTab = outTab[outTab[,"P-value"] < 0.05, ]
outTab = as.matrix(outTab[, 1:(ncol(outTab)-3)])
outTab = rbind(id = colnames(outTab), outTab)

# Save filtered CIBERSORT results
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)


# --- Correlation Analysis between Genes and Immune Cells ---

# Load required libraries
library(limma)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(linkET)

expFile = "normalized_data_matrix.txt"         # Normalized expression matrix
geneFile = "candidate_biomarkers.txt"          # Selected biomarker list
immFile = "CIBERSORT-Results.txt"              # CIBERSORT output file

# setwd("./")

# Read and format the expression data
rt = read.table(expFile, header=T, sep="\t", check.names=F)
rt = as.matrix(rt)
rownames(rt) = rt[,1]
exp = rt[, 2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data = avereps(data)

# Extract expression profiles for candidate biomarkers
geneRT = read.table(geneFile, header=F, sep="\t", check.names=F)
data = data[intersect(as.vector(geneRT[,1]), rownames(data)), ]

# Extract grouping info and retain only 'Treat' group samples
# Fixed the regex capture group backreference to \\2
group = gsub("(.*)\\_(.*?)", "\\2", colnames(data))
data = data[, group == "Treat", drop=F]
data = t(data)

# Read and format the immune infiltration results
immune = read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)

# Match samples between expression data and immune data
sameSample = intersect(row.names(data), row.names(immune))
data = data[sameSample, , drop=F]
immune = immune[sameSample, , drop=F]

# Remove immune cell types with zero variance
immune = immune[, apply(immune, 2, sd) > 0]

# Perform Spearman correlation analysis between genes and immune cells
geneCor = data.frame()
for(cell in colnames(immune)){
  for(gene in colnames(data)){
    x = as.numeric(immune[, cell])
    y = as.numeric(data[, gene])
    corT = cor.test(x, y, method="spearman")
    cor_val = corT$estimate
    p_val = corT$p.value
    geneCor = rbind(geneCor, cbind(spec=gene, env=cell, r=cor_val, p=p_val))
  }
}

# Format correlation results for visualization
geneCor$r = as.numeric(geneCor$r)
geneCor$p = as.numeric(geneCor$p)
geneCor$pd = ifelse(geneCor$p < 0.05, ifelse(geneCor$r > 0, "Postive", "Negative"), "Not")
geneCor$absR = abs(geneCor$r)

# Bin the correlation coefficients for aesthetic mapping
geneCor = geneCor %>% mutate(rd = cut(absR, breaks = c(-Inf, 0.2, 0.4, 0.6, Inf),
                                    labels = c("< 0.2", "0.2 - 0.4", "0.4 - 0.6", ">= 0.6")))

# Generate the integrated correlation plot using linkET
qcorPlot = qcorrplot(correlate(immune, method="spearman"), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), 
              data = geneCor, 
              curvature = nice_curvature()) +
  # Customize color gradients and legend titles
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu"))) +
  scale_size_manual(values = c(0.5, 1.5, 2, 3)) +
  scale_colour_manual(values = c("#1B9E77", "#CCCCCC99", "#D95F02")) +
  guides(size = guide_legend(title = "abs(Cor)",
                               override.aes = list(colour = "grey35"), 
                               order = 2),
         colour = guide_legend(title = "Significance", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Cell-cell cor", order = 3))

# Export the final visualization to a PDF file
pdf(file="linkET_correlation.pdf", width=9, height=7)
print(qcorPlot)
dev.off()