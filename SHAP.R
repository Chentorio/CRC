library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(kernelshap)
library(pROC)
library(shapviz)
library(xgboost)
library(klaR)

# Define input files
expFile = "normalized_data_matrix.txt"            # Normalized expression matrix
geneFile = "candidate_biomarkers.txt"             # Selected candidate biomarkers
# setwd("./")                                     # Set your working directory here

# Read the expression data
data = read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

# Read candidate biomarker list and filter the data
geneRT = read.table(geneFile, header=F, sep="\t", check.names=F)
validGenes = intersect(as.vector(geneRT[,1]), row.names(data))
data = data[validGenes,]

# Replace hyphens with underscores to avoid formula syntax errors
row.names(data) = gsub("-", "_", row.names(data))

# Transpose data and extract grouping information from sample names
data = t(data)
group = gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(data))
data = as.data.frame(data)
data$Type = group

# Partition data into training (70%) and testing (30%) sets
inTrain = createDataPartition(y=data$Type, p=0.7, list=F)
train = data[inTrain,]
test = data[-inTrain,]

# Prepare test set labels for evaluation
yTestClass = test$Type
yTest = ifelse(yTestClass == "Control", 0, 1)
test = test[, -ncol(test)]

# Set cross-validation parameters
control = trainControl(method="repeatedcv", number=5, savePredictions=TRUE)
methodRT = read.table("refer.methodLists.txt", header=T, sep="\t", check.names=F)

# Iterate through machine learning methods and evaluate via ROC
modelList = list()
AUC = c()
ROCcolor = rainbow(nrow(methodRT))

for(i in 1:nrow(methodRT)){
    name = methodRT[i, "Name"]
    method = methodRT[i, "Method"]
    
    # Train the model based on the specified method
    if(name == "SVM"){
        # SVM requires probability modeling to be explicitly enabled
        model = train(Type ~ ., data = train, method=method, prob.model=TRUE, trControl = control)
    } else {
        model = train(Type ~ ., data = train, method=method, trControl = control)
    }
    
    # Predict probabilities on the test set
    pred = predict(model, newdata=test, type="prob")
    
    # Calculate and store ROC information
    roc_obj = roc(yTest, as.numeric(pred[,2]))
    AUC = c(AUC, paste0(name, ": ", sprintf("%.03f", roc_obj$auc)))
    modelList[[method]] = as.numeric(roc_obj$auc)
    
    # Plot combined ROC curves
    if(i == 1){
        pdf(file="ROC.pdf", width=5.5, height=5)
        plot(roc_obj, print.auc=F, legacy.axes=T, main="", col=ROCcolor[i], lwd=3)
    } else {
        plot(roc_obj, print.auc=F, legacy.axes=T, main="", col=ROCcolor[i], lwd=3, add=T)
    }
}
legend("bottomright", AUC, col=ROCcolor, lwd=3, bty = "n", cex=0.9)
dev.off()

# Identify and select the model with the highest AUC
aucValue = unlist(modelList)
bestMethod = names(which(aucValue == max(aucValue)))

# Recast target variable to binary for SHAP calculation
train$Type = ifelse(train$Type == "Control", 0, 1)
# Use the first element in case of a tie in AUC
finalModel = train(Type ~ ., data = train, method = bestMethod[1], trControl=control)

# Calculate SHAP values using Permutation SHAP
# Note: kernelshap can be used for more complex models but takes longer
fit = permshap(finalModel, train[, -ncol(train)])
shp = shapviz(fit, X_pred = train[, -ncol(train)], X = train[, -ncol(train)], interactions=T)

# Sort features by importance based on mean absolute SHAP values
important = sort(colMeans(abs(shp$S)), decreasing=T)
showVars = names(important)

# Generate and save SHAP importance bar plot
pdf(file="barplot.pdf", width=6, height=6)
sv_importance(shp, kind="bar", show_numbers=TRUE) + theme_bw()
dev.off()

# Generate and save SHAP summary bee-swarm plot
pdf(file="bee.pdf", width=7, height=6)
sv_importance(shp, kind = "bee", show_numbers=TRUE) + theme_bw()
dev.off()

# Generate and save SHAP dependence plots for top variables
pdf(file="dependence.pdf", width=9, height=6)
sv_dependence(shp, v = showVars) + theme_bw()
dev.off()

# Generate and save SHAP waterfall plot for the first sample
pdf(file="waterfall.pdf", width=7, height=5)
sv_waterfall(shp, row_id = 1)
dev.off()

# Generate and save SHAP force plot for the first sample
pdf(file="force.pdf", width=9, height=5)
sv_force(shp, row_id = 1)
dev.off()