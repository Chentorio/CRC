set.seed(12345)
library(limma)
library(e1071)

# Define input file path
expFile = "diffGeneExp.txt"

# Set working directory (Relative path is recommended for portability)
# setwd("./")

# Load external SVM-RFE functions
source("msvmRFE.R")

# Read input expression data
rt = read.table(expFile, header=TRUE, sep="\t", check.names=FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[,1]
exp = rt[, 2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)

# Average duplicated gene expressions and transpose the matrix
data = avereps(data)
data = as.data.frame(t(data))

# Extract sample grouping information
group = gsub("(.*)_(.*)", "\\2", row.names(data))
data = cbind(group, data)
data$group = factor(data$group, levels=c("Control", "Treat"))

# Execute SVM-RFE algorithm for feature selection
# This process includes k-fold cross-validation
svmRFE(data, k=10, halve.above=50)
nfold = 10
geneNum = nrow(data)
folds = rep(1:nfold, len=geneNum)[sample(geneNum)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, data, k=10, halve.above=50)

# Extract and save ranked features based on importance
top.features = WriteFeatures(results, data, save=FALSE)
write.table(top.features, file="feature_svm.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Perform cross-validation to determine the optimal number of features
# Limit the sweep to the top 30 features if the total count is high
num = ifelse((ncol(data)-1) > 30, 30, (ncol(data)-1))
featsweep = lapply(1:num, FeatSweep.wrap, results, data)

# Calculate error rates from cross-validation
no.info = min(prop.table(table(data[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

# Plot and save error rate visualization
pdf(file="errors.pdf", width=5, height=5)
PlotErrors(errors, no.info=no.info)
dev.off()

# Plot and save accuracy visualization
pdf(file="accuracy.pdf", width=5, height=5)
Plotaccuracy(1-errors, no.info=no.info)
dev.off()

# Extract and save the final set of characteristic genes
featureGenes = top.features[1:which.min(errors), 1, drop=FALSE]
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)