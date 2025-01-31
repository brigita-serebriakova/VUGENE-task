#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("limma")

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("edgeR")

#install.packages("poorman")
#install.packages("ggplot2")

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("biomaRt")

#install.packages("ggrepel")

library(limma)
library(edgeR)
library(poorman)
library(ggplot2)
library(biomaRt)
library(ggrepel)
library(dplyr)

setwd('/Users/brigitaserebriakova/Desktop/VUGENE-task')

counts_matrix <- read.csv('rnaseq_default_counts.csv', header = T, row.names = 1)
covariates <- read.csv('clean_sample.csv', header = T)

all(colnames(counts_matrix) == covariates$SID)

# Creating a design matrix

covariates$Condition <- factor(covariates$Condition, levels = c("Control", "AlcoholUseDisorder"))
covariates$Batch <- factor(covariates$Batch)
covariates$Sex <- factor(covariates$Sex)
covariates$Liver <- factor(covariates$Liver, levels = c("Normal", "Steatosis", "Cirrhosis", "Congestion"))
covariates$LeftRightBrain <- factor(covariates$LeftRightBrain)
covariates$Smoking <- factor(covariates$Smoking, levels = c("Current", "Never", "Past", "Unknown"))


design <- model.matrix(~ Condition + Batch + Sex + PMI + RIN + Liver + BrainPH + BrainWeight + LeftRightBrain + Smoking, data = covariates) # look into it
#colnames(design) <- c("Control", "AUD")

### Normalization

dge <- DGEList(counts = counts_matrix)

# Filtering out genes that are not expressed
keep <- filterByExpr(dge, design)
dge <- dge[keep,]
nrow(dge$counts)

# Filter out genes that fail to show at least 1 count-per-million (cpm) reads in at least 12 samples
isexpr <- rowSums(cpm(dge) > 1) >= 12
dge <- dge[isexpr,,keep.lib.sizes=FALSE]
nrow(dge$counts)

# TMM normalization
dge <- calcNormFactors(dge)


### Quality control

low_RIN <- c("SRR15466722", "SRR15466724", "SRR15466742", "SRR15466727", "SRR15466740")

# All genes by batch
#png("by-batch.png")
plotMDS(dge, pch = 19, col = ifelse(covariates$Batch == 1, "blue", "red"))
text(x = plotMDS(dge, plot = FALSE)$x[covariates$SID %in% low_RIN],
     y = plotMDS(dge, plot = FALSE)$y[covariates$SID %in% low_RIN],
     labels = covariates$SID[covariates$SID %in% low_RIN],
     pos = 4, cex = 0.8)
#dev.off()

#png("by-batch-corrected.png")
dge2 <- removeBatchEffect(dge)
plotMDS(dge2, pch = 19, col = ifelse(covariates$Batch == 1, "blue", "red"))
#dev.off()

# All genes by Condition
#png("by-condition.png")
plotMDS(dge, pch = 19, col = ifelse(covariates$Condition == "Control", "blue", "red"))
#dev.off()

# Top 100 genes by batch
#png("by-batch-top100.png")
plotMDS(dge, top = 100, pch = 19, col = ifelse(covariates$Batch == 1, "blue", "red"))
#dev.off()

# Top 100 genes by Condition
#png("by-condition-top100.png")
plotMDS(dge, top = 100, pch = 19, col = ifelse(covariates$Condition == "Control", "blue", "red"))
#dev.off()

# Top 100 genes by Sex
#png("by-sex-top100.png")
plotMDS(dge, top = 100, pch = 19, col = ifelse(covariates$Sex == "Male", "blue", "red"))
#dev.off()


# voom

v <- voom(dge, design)#, normalize.method = "quantile") #check page 127
fit <- lmFit(v, design) #look into correlation parameter?
contrast_matrix <- makeContrasts(AUDvsControl = ConditionAlcoholUseDisorder - 0, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2<- eBayes(fit2)
summary(decideTests(fit2)) # what does it show?
results <- topTable(fit2, number = Inf)


# marking upregulated, downregulated and non-significant genes
results$direction <- with(results, poorman::case_when(
  logFC > 1 & adj.P.Val < 0.05 ~ "up",
  logFC < -1 & adj.P.Val < 0.05 ~ "down",
  TRUE ~ "n.s."
))

results$direction = factor(results$direction, levels = c("down", "up", "n.s."))

volcanoplot(fit2)

# finding gene symbols for top 10 genes
top10 <- topTable(fit2)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ens <- rownames(top10)
gene_symbols <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = ens, mart = ensembl)

results$ENS <- rownames(results)

results <- results %>%
  left_join(gene_symbols, by = c("ENS" = "ensembl_gene_id"))

labels <- results %>% filter(!is.na(hgnc_symbol))


# Volcano plot with top 10 genes
volc_plot2 <- ggplot() +
  geom_point(data = results,
            aes(x = logFC, y = -log10(P.Value), color = direction),
            alpha = 0.5) +
  scale_color_manual(values = c("up" = "#E41A1C",
                                "down" = "#377EB8",
                                "n.s." = "lightgrey"),
                     guide = "none") +
  geom_text_repel(data = labels,
                  aes(x = logFC, y = -log10(P.Value), label = hgnc_symbol),
                  size = 3,
                  box.padding = 0.1) +
  labs(
    x = "Fold change (log2)",
    y = "-log10(p-value)"
  ) +
  theme(panel.grid = element_blank()) +
  geom_hline(yintercept = max(-log10(results$P.Value)) - 3.85, linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")
volc_plot2

#png("volcano-plot.png")
#print(volc_plot2)
#dev.off()




volc_plot <- ggplot() +
  geom_point(data = results,
            aes(x = logFC, y = -log10(P.Value), color = direction),
            alpha = 0.5) +
  scale_color_manual(values = c("up" = "#E41A1C",
                                "down" = "#377EB8",
                                "n.s." = "lightgrey"),
                     guide = "none") +
  labs(
    x = "Fold change (log2)",
    y = "-log10(p-value)"
  ) +
  theme(panel.grid = element_blank()) +
  geom_hline(yintercept = max(-log10(results$P.Value)) - 3.85, linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")
volc_plot

### older analysis
counts_matrix <- read.csv('rnaseq_default_counts.csv', header = T, row.names = 1)
covariates <- read.csv('clean_sample.csv', header = T)

covariates$Condition <- factor(covariates$Condition, levels = c("Control", "AlcoholUseDisorder"))
covariates$Batch <- factor(covariates$Batch)

design <- model.matrix(~ Batch + Condition, data = covariates) # look into it

dge <- DGEList(counts = counts_matrix)

keep <- filterByExpr(dge, design)
dge <- dge[keep,]
nrow(dge$counts)

isexpr <- rowSums(cpm(dge) > 1) >= 12
dge <- dge[isexpr,,keep.lib.sizes=FALSE]
nrow(dge$counts)

dge <- calcNormFactors(dge)

v <- voom(dge, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)
summary(decideTests(fit2)) # what does it show?
topTable(fit)

#png("first-volcano-plot.png")
#volcanoplot(fit)
#dev.off()




