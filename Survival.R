library(limma)
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(dplyr)

# Global Settings
target_gene <- "PSAT1"
exp_file <- "diffGeneExp.txt"      # Normalized expression matrix
time_file <- "time.txt"           # Survival data (futime, fustat)
cli_file <- "clinical.txt"         # Clinical metadata (age, stage, etc.)
diff_res_file <- "diff.txt"        # Differential analysis results (for logFC check)

# Set working directory
# setwd("./") 

# Load expression data
rt = read.table(exp_file, header=T, sep="\t", check.names=F)
rt = as.matrix(rt)
rownames(rt) = rt[,1]
exp = rt[, 2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data = avereps(data)
data = log2(data + 1) # Normalize if not already logged

# Filter Tumor samples (TCGA logic: sample code 01-09 is tumor)
group = sapply(strsplit(colnames(data), "\\-"), "[", 4)
group = sapply(strsplit(group, ""), "[", 1)
data = data[, group == 0]

# Trim barcodes to Patient ID (e.g., TCGA-XX-XXXX)
colnames(data) = gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data))
data = t(data)
data = avereps(data)

# Load and merge survival data
cli_surv = read.table(time_file, header=T, sep="\t", check.names=F, row.names=1)
same_samples = intersect(row.names(data), row.names(cli_surv))
data = data[same_samples, ]
cli_surv = cli_surv[same_samples, ]
merged_data = cbind(cli_surv, data)


# Convert days to years
merged_data$futime = merged_data$futime / 365
diff_res = read.table(diff_res_file, header=T, sep="\t", check.names=F, row.names=1)

surv_results = data.frame()
# Loop through genes (skipping futime and fustat columns)
for(gene in colnames(merged_data)[3:ncol(merged_data)]){
    if(sd(merged_data[,gene]) < 0.1){ next }
    
    # Univariate Cox Analysis
    cox_fit = coxph(Surv(futime, fustat) ~ merged_data[,gene], data = merged_data)
    cox_sum = summary(cox_fit)
    log_hr = log2(cox_sum$conf.int[,"exp(coef)"])
    cox_p = cox_sum$coefficients[,"Pr(>|z|)"]
    
    # KM Analysis (Median split)
    group_km = ifelse(merged_data[,gene] > median(merged_data[,gene]), "high", "low")
    km_diff = survdiff(Surv(futime, fustat) ~ group_km, data = merged_data)
    km_p = 1 - pchisq(km_diff$chisq, df=1)
    
    # Selection criteria: P < 0.05 and consistency with differential logFC
    if((cox_p < 0.05) & (km_p < 0.05)){
        if(gene %in% rownames(diff_res) && (log_hr * diff_res[gene, "logFC"] > 0)){
            surv_results = rbind(surv_results, cbind(id=gene, HR=cox_sum$conf.int[,"exp(coef)"], 
                                                     coxP=cox_p, kmP=km_p))
        }
    }
}
write.table(surv_results, file="prognostic_genes.txt", sep="\t", row.names=F, quote=F)


clinical_meta = read.table(cli_file, header=T, sep="\t", check.names=F, row.names=1)
same_cli = intersect(row.names(merged_data), row.names(clinical_meta))
exp_sub = merged_data[same_cli, ]
cli_sub = clinical_meta[same_cli, ]

indep_results = data.frame()
for(i in as.vector(surv_results$id)){
    if(!(i %in% colnames(exp_sub))){ next }
    rt_multi = cbind(exp_sub[,1:2], cli_sub, gene=exp_sub[,i])
    multi_cox = coxph(Surv(futime, fustat) ~ ., data = rt_multi)
    multi_sum = summary(multi_cox)
    
    # Filter for genes that remain significant independent of clinical factors
    if(multi_sum$coefficients["gene","Pr(>|z|)"] < 0.05){
        indep_results = rbind(indep_results, cbind(id=i, HR=multi_sum$conf.int["gene","exp(coef)"], 
                                                   pvalue=multi_sum$coefficients["gene","Pr(>|z|)"]))
    }
}
write.table(indep_results, file="independent_prognostic_genes.txt", sep="\t", row.names=F, quote=F)


# Use original mixed data (Tumor + Normal) for boxplot
rt_box = read.table(exp_file, header=T, sep="\t", check.names=F, row.names=1)
target_data = t(rt_box[target_gene, , drop=F])
group_box = sapply(strsplit(rownames(target_data),"\\-"), "[", 4)
group_box = ifelse(sapply(strsplit(group_box,""), "[", 1) == "0", "Tumor", "Normal")

df_box = data.frame(Expression = log2(as.numeric(target_data) + 1), Type = group_box)

pdf(file=paste0(target_gene, "_diff_boxplot.pdf"), width=5, height=4.5)
ggboxplot(df_box, x="Type", y="Expression", color="Type", palette = c("blue","red"), add = "jitter") +
    stat_compare_means(comparisons=list(c("Normal", "Tumor")), label = "p.signif") +
    ylab(paste(target_gene, "expression"))
dev.off()

# Find optimal cutoff for PSAT1 in tumor samples
res_cut = surv_cutpoint(merged_data, time = "futime", event = "fustat", variables = target_gene)
res_cat = surv_categorize(res_cut)
merged_data$Group = res_cat[, target_gene]

fit = survfit(Surv(futime, fustat) ~ Group, data = merged_data)
p_val = surv_pvalue(fit)$p.txt

surv_plot = ggsurvplot(fit, data=merged_data, pval=p_val, legend.title=target_gene,
                       legend.labs=c("High", "Low"), xlab="Time(years)", 
                       palette=c("#E28077FF", "#7BA0D8FF"), risk.table=T)

pdf(file=paste0(target_gene, "_optimal_survival.pdf"), width=5.5, height=5, onefile=F)
print(surv_plot)
dev.off()