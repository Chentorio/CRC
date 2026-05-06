set.seed(12345)
library(glmnet)

# Define input file
inputFile = "diffGeneExp.txt"

# Set working directory to the current path
# setwd("./")

# Read input file and transpose the matrix
rt = read.table(inputFile, header=TRUE, sep="\t", check.names=FALSE, row.names=1)
rt = t(rt)

# Prepare data for modeling
x = as.matrix(rt)
# Extract group information and convert to factor
y_char = gsub("(.*)_(.*)", "\\2", row.names(rt))
y = as.factor(y_char)

# Build LASSO regression model
fit = glmnet(x, y, family = "binomial", alpha=1)

# Plot Lasso regression coefficients and save to PDF
pdf(file="lasso.pdf", width=6, height=5.5)
plot(fit)
dev.off()

# Perform cross-validation and save the plot
cvfit = cv.glmnet(x, y, family="binomial", alpha=1, type.measure='deviance', nfolds = 10)
pdf(file="cvfit.pdf", width=6, height=5.5)
plot(cvfit)
dev.off()

# Extract and save selected genes based on optimal lambda
coef = coef(fit, s=cvfit$lambda.min)
index = which(coef != 0)
lassoGene = row.names(coef)[index]

# Remove the intercept from the gene list
lassoGene = lassoGene[-1]

# Output the results to a text file
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)