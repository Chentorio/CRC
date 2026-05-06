suppressPackageStartupMessages({
  library(limma)
  library(Seurat)
  library(dplyr)
  library(magrittr)
  library(celldex)
  library(SingleR)
  library(monocle)
  library(clustree)
  library(harmony)
  library(assertthat)
  library(SCpubr)
  library(patchwork)
})

# Verify that key packages are loaded from the custom library
check_pkg <- function(pkg) {
  p <- find.package(pkg)
  message(pkg, " -> ", p)
  if (!startsWith(p, my_lib)) {
    warning(pkg, " is not loaded from my_lib. Path=", p)
  }
}
check_pkg("Seurat")
check_pkg("ggplot2")

## ===================== Parameters =====================
logFCfilter <- 1
adjPvalFilter <- 0.05
workDir <- "/dssg/home/acct-medlf/medlf05/cmm/knockdown"
setwd(workDir)

## ===================== 01. Read 10X Genomics Data =====================
# Identify directories containing 10X matrices
all_dirs <- list.dirs(workDir, full.names = TRUE, recursive = FALSE)

is_10x_dir <- function(d) {
  file.exists(file.path(d, "filtered_feature_bc_matrix")) ||
    file.exists(file.path(d, "filtered_gene_bc_matrices")) ||
    file.exists(file.path(d, "matrix.mtx")) ||
    file.exists(file.path(d, "matrix.mtx.gz"))
}

dirs_sample <- all_dirs[vapply(all_dirs, is_10x_dir, logical(1))]
if (length(dirs_sample) == 0) {
  stop("No 10X directories found. Check folder structure at: ", workDir)
}

names(dirs_sample) <- basename(dirs_sample)

# Load data and create Seurat object
counts <- Read10X(data.dir = dirs_sample)
pbmc <- CreateSeuratObject(counts, min.cells = 3, min.features = 100)

## ===================== Quality Control (QC) =====================
# Calculate mitochondrial gene percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pdf("01.featureViolin.pdf", width = 10, height = 6.5)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Subset cells (Suggested: nFeature_RNA > 200 for better cell quality)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 15)

pdf("01.featureCor.pdf", width = 13, height = 7)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 1.5)
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 1.5)  
print(plot1 + plot2)
dev.off()

## ===================== Normalization & Variable Features =====================
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 1500)

top10 <- head(VariableFeatures(pbmc), 10)

pdf("01.featureVar.pdf", width = 10, height = 6)
p1 <- VariableFeaturePlot(pbmc)
p2 <- LabelPoints(plot = p1, points = top10, repel = TRUE)
print(p1 + p2)
dev.off()

## ===================== 02. PCA & Harmony Integration =====================
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))
pbmc <- RunPCA(pbmc, npcs = 20, features = VariableFeatures(pbmc))

# Correct batch effects across samples using Harmony
pbmc <- RunHarmony(pbmc, group.by.vars = "orig.ident")

# PCA Visualizations
pdf("02.pcaGene.pdf", width = 10, height = 8)
VizDimLoadings(pbmc, dims = 1:4, reduction = "pca", nfeatures = 20)
dev.off()

pdf("02.PCA.pdf", width = 7.5, height = 5)
DimPlot(pbmc, reduction = "pca")
dev.off()

pdf("02.pcaHeatmap.pdf", width = 10, height = 8)
DimHeatmap(pbmc, dims = 1:4, cells = 500, balanced = TRUE, nfeatures = 30, ncol = 2)
dev.off()

# JackStraw for dimensionality assessment
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
pdf("02.pcaJackStraw.pdf", width = 8, height = 6)
JackStrawPlot(pbmc, dims = 1:20)
dev.off()

## ===================== 03. Clustering & Marker Detection =====================
pcSelect <- 20
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:pcSelect)
pbmc <- FindClusters(pbmc, resolution = 0.6)
pbmc <- RunTSNE(pbmc, reduction = "harmony", dims = 1:pcSelect)

pdf("03.cluster.pdf", width = 7, height = 6)
print(TSNEPlot(pbmc, pt.size = 2, label = TRUE))
dev.off()

# Find markers for each cluster
pbmc.markers <- FindAllMarkers(
  pbmc,
  only.pos = FALSE,
  min.pct = 0.25,
  logfc.threshold = logFCfilter
)

sig.markers <- pbmc.markers[
  (abs(as.numeric(pbmc.markers$avg_log2FC)) > logFCfilter &
     as.numeric(pbmc.markers$p_val_adj) < adjPvalFilter),
]

write.table(sig.markers, file = "03.clusterMarkers.txt", sep = "\t", row.names = FALSE, quote = FALSE)

## ===================== 04. SingleR Cell Type Annotation =====================
# Extract data for SingleR (Seurat v5 uses layer syntax)
pbmc_for_SingleR <- GetAssayData(pbmc, layer = "data")
clusters <- pbmc@meta.data$seurat_clusters

# Load reference datasets (Ensure these files are in the workDir)
ref1 <- get(load("ref_Human_all.RData"))
ref2 <- get(load("ref_Hematopoietic.RData"))
ref3 <- get(load("DatabaseImmuneCellExpressionData.Rdata"))
ref4 <- get(load("BlueprintEncode_bpe.se_human.RData"))
ref5 <- get(load("HumanPrimaryCellAtlas_hpca.se_human.RData"))
ref6 <- get(load("MonacoImmuneData.Rdata"))
ref7 <- get(load("NovershternHematopoieticData.Rdata"))

# Predict cell types using multiple references
singler <- SingleR(
  test = pbmc_for_SingleR,
  ref = list(ref1, ref2, ref3, ref4, ref5, ref6, ref7),
  labels = list(ref1$label.main, ref2$label.main, ref3$label.main, ref4$label.main,
                ref5$label.main, ref6$label.main, ref7$label.main),
  clusters = clusters
)

# Clean up cell type labels
singler$labels <- gsub("_|-", " ", singler$labels)
singler$labels <- gsub("T cells, CD4\\+", "CD4+ T cells", singler$labels)
singler$labels <- gsub("T cells, CD8\\+", "CD8+ T cells", singler$labels)

# Export annotation results
clusterAnn <- as.data.frame(singler)
clusterAnn <- cbind(id = row.names(clusterAnn), clusterAnn)
clusterAnn <- clusterAnn[, c("id", "labels")]
write.table(clusterAnn, file = "04.clusterAnn.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# Rename identities in Seurat object
newLabels <- singler$labels
names(newLabels) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, newLabels)

pdf("04.cellAnn.pdf", width = 7.5, height = 6)
print(TSNEPlot(pbmc, pt.size = 2, label = TRUE))
dev.off()

## ===================== Grouping: Control vs Treat =====================
# Assign groups based on sample name prefix
pbmc$Type <- ifelse(grepl("^treat", pbmc$orig.ident, ignore.case = TRUE), "Treat", "Control")

pdf("04.group.cellAnn.pdf", width = 11, height = 6)
print(TSNEPlot(pbmc, pt.size = 1, label = TRUE, split.by = "Type"))
dev.off()

## ===================== Treat vs Control DE per Cell Type =====================
cellAnn_vec <- clusterAnn[match(pbmc$seurat_clusters, clusterAnn[, 1]), 2]
pbmc$group <- paste0(pbmc$Type, "_", cellAnn_vec)

for (cellName in unique(cellAnn_vec)) {
  controlName <- paste0("Control_", cellName)
  treatName   <- paste0("Treat_", cellName)
  
  # Minimum requirement: at least 5 cells per group
  if ((sum(pbmc$group == controlName) > 5) && (sum(pbmc$group == treatName) > 5)) {
    de <- FindMarkers(
      pbmc,
      ident.1 = treatName,
      ident.2 = controlName,
      group.by = "group",
      logfc.threshold = 0.1
    )
    
    sig.de <- de[(abs(de$avg_log2FC) > logFCfilter & de$p_val_adj < adjPvalFilter), , drop = FALSE]
    sig.de <- cbind(Gene = row.names(sig.de), sig.de)
    
    # Export DE results with cleaned filenames
    file_name <- paste0("05.", gsub("[/\\\\: ]+", "_", cellName), ".Treat_vs_Control.diffGene.txt")
    write.table(sig.de, file = file_name, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

## ===================== Save Results =====================
save(pbmc, sig.markers, file = "Seurat_Final.Rdata")