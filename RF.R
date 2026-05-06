set.seed(12345)
library(randomForest)
library(ggplot2)

# Define the input file name
inputFile = "diffGeneExp.txt"

# Set working directory to the current path
# setwd("./")

# Load and prepare data
data = read.table(inputFile, header=TRUE, sep="\t", check.names=FALSE, row.names=1)
data = t(data)

# Extract group information from row names
group = gsub("(.*)_(.*)", "\\2", row.names(data))

# Build initial Random Forest model with 500 trees
rf = randomForest(as.factor(group) ~ ., data=data, ntree=500)

# Save the error rate plot to a PDF
pdf(file="forest.pdf", width=6, height=6)
plot(rf, main="Random forest", lwd=2)
dev.off()

# Identify the tree number with the minimum OOB error
optionTrees = which.min(rf$err.rate[,1])
print(paste("Optimized number of trees:", optionTrees))

# Refit the Random Forest model using the optimized tree number
rf2 = randomForest(as.factor(group) ~ ., data=data, ntree=optionTrees)

# Retrieve feature importance scores
importance = importance(x=rf2)

# Sort and filter genes based on Mean Decrease Gini score
rfGenes_all = importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
# Filter genes with an importance score greater than 2
rfGenes = names(rfGenes_all[rfGenes_all > 2]) 

# Save the filtered gene list
write.table(rfGenes, file="RF.gene.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

# Prepare data for the importance bubble plot
importance_sorted = importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE), , drop=FALSE]
rt = cbind(row.names(importance_sorted), as.data.frame(importance_sorted))
colnames(rt) = c("ID", "Importance")

# Determine the number of genes to display for the plot
showNum = 30
if(nrow(rt) < showNum){
    showNum = nrow(rt)
}
rt = rt[1:showNum,]

# Format data for ggplot2
rt$ID = factor(rt$ID, levels=rev(rt$ID))
rt$Importance = as.numeric(rt$Importance)

# Create the importance visualization
p = ggplot(rt, aes(Importance, ID)) + 
    geom_segment(aes(x=0, xend=Importance, y=ID, yend=ID), color="grey", cex=1.5) + 
    geom_point(aes(color=Importance), cex=3.5) +
    scale_colour_gradient(low="#008280FF", high="#BB0021FF") + 
    labs(x="Importance", y="") +
    theme_bw() +
    theme(axis.text.x=element_text(color="black", size=10),
          axis.text.y=element_text(color="black", size=10))

# Export the plot to a PDF
pdf(file="geneImportance.pdf", width=5, height=4.5)
print(p)
dev.off()