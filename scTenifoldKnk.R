suppressPackageStartupMessages({
  library(scTenifoldKnk)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(igraph)
  library(ggrepel)
})

set.seed(123)

## ===== User settings =====
target_gene <- "PSAT1"   # Target gene for virtual knockout
work_dir <- "/"
setwd(work_dir)

## ===== Helper: GetAssayData compatible for Seurat v4/v5 =====
GetCountsMatrix <- function(seurat_obj, assay = NULL) {
  if (!is.null(assay)) DefaultAssay(seurat_obj) <- assay
  # Compatibility check for Seurat v5 (layer) and v4 (slot)
  out <- tryCatch(GetAssayData(seurat_obj, layer = "counts"),
                  error = function(e) NULL)
  if (is.null(out)) {
    out <- GetAssayData(seurat_obj, slot = "counts")
  }
  return(out)
}

## ===== Load Seurat object =====
load("Seurat.Rdata")

if (!exists("pbmc")) {
  stop("Object pbmc not found in Seurat.Rdata. Please verify the object name.")
}

## ===== Subset to Treat group for analysis =====
if (!("Type" %in% colnames(pbmc@meta.data))) {
  message("Meta column Type not found. Inferring group from orig.ident.")
  pbmc$Type <- ifelse(grepl("^treat", pbmc$orig.ident, ignore.case = TRUE), "Treat", "Control")
}

message("Current Group Distribution:")
print(table(pbmc$Type))

pbmc <- subset(pbmc, subset = Type == "Treat")

if (ncol(pbmc) < 50) {
  warning("Insufficient cell count in Treat group (n=", ncol(pbmc), "). Results may be unstable.")
}

## ===== Extract and validate count matrix =====
countMat <- GetCountsMatrix(pbmc)

# Ensure the target gene is present in the dataset
if (!(target_gene %in% rownames(countMat))) {
  stop("Target gene ", target_gene, " not found in matrix. Check gene symbols.")
}

## ===== Select features (HVGs + Target Gene) =====
# Note: 10,000 features may require significant RAM on HPC
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 10000)
hvgs <- VariableFeatures(pbmc)

genes_use <- unique(c(target_gene, hvgs))
genes_use <- intersect(genes_use, rownames(countMat)) 

data <- as.data.frame(as.matrix(countMat[genes_use, , drop = FALSE]))

## ===== Execute Virtual Knockout Pipeline =====
# nc_nNet: Number of networks; nc_nCells: Cells per network subsampling
result <- scTenifoldKnk(
  countMatrix = data,
  gKO = target_gene,
  qc = TRUE,
  qc_mtThreshold = 0.1,
  qc_minLSize = 1000,
  nc_nNet = 10,
  nc_nCells = 500
)

## ===== Export Differential Regulation Results =====
df <- result$diffRegulation
df <- df[df$gene != target_gene, , drop = FALSE]

write.table(df, file = "diffRegulation_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Filter for significant regulatory changes (p.adj < 0.05)
outTab <- df[df$p.adj < 0.05, , drop = FALSE]
write.table(outTab, file = "sigDiff.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## =========================== Barplot ===========================
# Visualize the top 20 genes with the highest regulatory disruption
top_genes <- head(df[order(-df$FC), ], 20)

p1 <- ggplot(top_genes, aes(x = reorder(gene, FC), y = FC)) +
  geom_bar(stat = "identity", fill = "#5A9BD4") +
  coord_flip() +
  labs(title = "Top 20 Differentially Regulated Genes", x = "Gene", y = "Regulatory FC") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

pdf(file = "barplot.pdf", width = 6, height = 5)
print(p1)
dev.off()

## =========================== Volcano Plot ===========================
df$log_p.adj <- -log10(df$p.adj)
df$significant <- ifelse(df$p.adj < 0.05, "Significant", "Not significant")
label_genes <- subset(df, p.adj < 0.05)

# Calculate a cap for the Y-axis to maintain visual clarity
y_upper <- quantile(df$log_p.adj, 0.999, na.rm = TRUE)

p2 <- ggplot(df, aes(x = Z, y = log_p.adj, color = significant)) +
  geom_point(alpha = 0.7, size = 1) +
  scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray50")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(data = label_genes, aes(label = gene), size = 3, max.overlaps = 50) +
  labs(title = "Knockout Impact Distribution", x = "Z-score", y = "-log10(p.adj)") +
  theme_classic() +
  coord_cartesian(ylim = c(0, y_upper)) +
  theme(legend.position = "none")

pdf(file = "vol.pdf", width = 6, height = 5)
print(p2)
dev.off()

## =========================== Proportion Pie Chart ===========================
sig_count <- table(df$significant)
sig_df <- as.data.frame(sig_count)
colnames(sig_df) <- c("category", "count")
sig_df$percentage <- paste0(round(sig_df$count / sum(sig_df$count) * 100, 1), "%")

p3 <- ggplot(sig_df, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = percentage), position = position_stack(vjust = 0.5), size = 4) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("Significant" = "red", "Not significant" = "lightgray")) +
  labs(title = "Regulatory Disruption Proportion", fill = "Significance") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

pdf(file = "pie.pdf", width = 6, height = 5)
print(p3)
dev.off()