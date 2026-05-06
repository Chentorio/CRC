library(dplyr)
library(ggplot2)
library(ggrepel)

logFCfilter = 0.585                   # logFC threshold
adj.P.Val.Filter = 0.05               # Adjusted p-value threshold
diffFile = "differential_expression_results.txt"  # Input file for all differential results
geneFile = "candidate_biomarkers.txt"             # File containing target gene list

# Set working directory to the current path
# setwd("./")

# Read the differential expression results
rt = read.table(diffFile, header=TRUE, sep="\t", check.names=FALSE)
row.names(rt) = rt[,1]

# Categorize genes into Up, Down, or Not significant
Sig = ifelse((rt$adj.P.Val < adj.P.Val.Filter) & (abs(rt$logFC) > logFCfilter), 
             ifelse(rt$logFC > logFCfilter, "Up", "Down"), "Not")

# Prepare data for plotting
rt = mutate(rt, Sig = Sig)

# Create the basic volcano plot
# Removed xlim to prevent data loss
p = ggplot(rt, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(col = Sig)) +
    scale_color_manual(values = c("green", "grey", "red")) +
    labs(title = " ") +
    theme_bw() + # Apply theme_bw here for a clean background
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))

# Label the candidate biomarkers on the plot
geneRT = read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)
sameGene = intersect(as.vector(geneRT[,1]), row.names(rt))
showData = rt[sameGene,]

# Add labels with repel to avoid overlapping
p1 = p + geom_label_repel(data = showData,
                         box.padding = 0.2, 
                         point.padding = 0.2, 
                         min.segment.length = 0.1,
                         size = 3, 
                         aes(label = id))

# Save the volcano plot to a PDF file
pdf(file = "volcano_plot.pdf", width = 5.25, height = 4.5)
print(p1)
dev.off()