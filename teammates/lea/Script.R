# Lea's analysis, exported from teammates/lea/Script.ipynb (R). Text cells appear as comments.
# Source: lhasch/System-Genomics-Project@8a133ad. Run from Lea's Data/ folder (see setwd() below).

# ====================================================================================================
# # Rejuvenating interventions and aging in the mouse cerebellum
# **Objective:** characterize the transcriptional effects of two rejuvenating interventions, acute dietary restriction (aDR) and young mouse plasma (YMP), and test how they relate to the aging signature.
#
# - **Reference:** Hahn et al. (2023)
# - **Data:** NCBI BioProjects PRJNA875066 (aging cohort, 12 runs) and PRJNA946653 (intervention cohort, 15 runs); run accessions in `samples.csv`.
# ====================================================================================================

# ====================================================================================================
# ## Workflow
# 1. Data import and sample exclusion
# 2. Sample overview: gene filtering, correlation clustering, PCA
# 3. Differential expression (DESeq2)
# 4. Gene set enrichment (GSEA)
# 5. Intervention effects vs the aging signature
# 6. Conclusions
# ====================================================================================================

# Install the R and Bioconductor packages used below (run once; skip if installed)
install.packages(c("tidyverse", "ggrepel", "BiocManager", "pbapply", "gplots", "msigdbr", "viridis", "languageserver", "httpgd"))
BiocManager::install(c("biomaRt", "sva", "DESeq2", "edgeR", "tximport", "fgsea", "apeglm"))

# Plain-text tables in Jupyter (no LaTeX / HTML rendering)
options(jupyter.rich_display = FALSE)

# ====================================================================================================
# ## 1 · Data import
# ====================================================================================================

# RSEM output, imported as gene × sample matrices (TPM and expected counts), with sample metadata and GENCODE M39 gene annotation.

# 1 · Load the data: TPM and expected counts (genes x samples), sample metadata, gene symbols.
# Adds a dataset label per sample and a cache helper for the slow steps.
# setwd() points to Lea's local Data folder: change it to run elsewhere.
setwd("C:/Users/hasch/OneDrive/Documents/BSSE/Systems Genomics/Data Analysis Project 25.11.2026/Data")
options(jupyter.rich_display = FALSE)
library(dplyr)
library(ggplot2)

# Expression data: genes x samples
expr   <- as.matrix(read.csv("rsem_tpm.csv", row.names = 1, check.names = FALSE))
counts <- as.matrix(read.csv("rsem_expected_counts.csv", row.names = 1, check.names = FALSE))
counts <- counts[rownames(expr), colnames(expr)]

# Sample metadata, same order as the matrix columns
meta <- read.csv("samples.csv")
meta <- meta[match(colnames(expr), meta$SRR), ]

# Label each sample with the dataset it comes from
meta$dataset <- ifelse(meta$cohort == "Aging", "PRJNA875066", "PRJNA946653")

# Set the order of the groups (controls / young first)
meta$group <- factor(meta$group, levels = c("aging 3mo", "aging 18mo", "aging 21mo",
                                            "AL", "aDR", "PBS", "YMP"))

# Gene metadata, same order as the matrix rows
genes <- read.csv("genes_M39.csv")
meta_genes <- data.frame(gene_id = rownames(expr),
                         symbol  = genes$symbol[match(rownames(expr), genes$gene_id)])

# Cache for the slow steps (DESeq2 fits, apeglm shrinkage, GSEA, MSigDB download)
# First run: each result is computed and saved to Data/cache/<name>.rds. Later runs: loaded from there.
# Delete the Data/cache folder (or set use_cache <- FALSE) after changing samples, filters or designs.
use_cache <- TRUE
dir.create("cache", showWarnings = FALSE)
cached <- function(name, expr) {
  f <- file.path("cache", paste0(name, ".rds"))
  if (use_cache && file.exists(f)) return(readRDS(f))
  x <- expr          # only computed when it is not loaded from the cache
  saveRDS(x, f)
  x
}

table(meta$dataset, meta$group)

# ====================================================================================================
# ### Sample exclusion
# **SRR21355389** (aging, 21 months, female) is excluded from all analyses: its haemoglobin genes (*Hbb-bs, Hba-a1, Hba-a2, Hbb-bt*) are about 15× higher than in any other sample, suggesting blood contamination.
# ====================================================================================================

# 1 · Sample exclusion: haemoglobin TPM per sample, then SRR21355389 is removed from all tables
# Hemoglobin genes (TPM) as a marker of blood contamination
hb_genes <- meta_genes$gene_id[meta_genes$symbol %in% c("Hbb-bs", "Hba-a1", "Hba-a2", "Hbb-bt")]
hb <- round(t(expr[hb_genes, ]))
colnames(hb) <- meta_genes$symbol[match(hb_genes, meta_genes$gene_id)]
rownames(hb) <- paste(meta$SRR, meta$group)
hb

# Exclude SRR21355389 (aging 21 mo, female) from all analyses
excluded <- "SRR21355389"
keep   <- !colnames(expr) %in% excluded
expr   <- expr[, keep]
counts <- counts[, keep]
meta   <- meta[keep, ]

table(meta$dataset, meta$group)

# Genes x samples after the exclusion
dim(expr)

# ====================================================================================================
# ## 2 · Sample overview
# **Expressed genes** (used in this section): detected in at least half of the samples, or mean TPM ≥ 1.
# ====================================================================================================

# 2 · Mean TPM per gene, raw and log scale: most genes are barely expressed
avg_expr <- rowMeans(expr)

layout(matrix(1:2, nrow=1))
hist(avg_expr)
hist(log10(avg_expr + 1))


library(ggplot2)

ggplot(data.frame(avg_expr), aes(x=avg_expr)) +
  geom_histogram(bins = 50) +
  scale_x_continuous(breaks = c(0,1,10,100,1000,10000,20000), trans="log1p", expand=c(0,0)) +
  scale_y_continuous(breaks = c(0,1), expand=c(0,0), trans="log1p") +
  theme_minimal()

# 2 · Expressed-gene filter for the sample overview: TPM > 0 in >= half the samples, or mean TPM >= 1.
# Genes are flagged (meta_genes$expressed), not removed, so DESeq2 later starts from all genes.
# Number of samples in which each gene is detected
num_det <- rowSums(expr > 0)
hist(num_det)

# Using the threshold of 1) being detected in at least half of the samples, or 2) average TPM>=1, to define whether a gene is expressed or not.
expressed <- rowMeans(expr > 0) >= 0.5 | rowMeans(expr) >= 1

## Removing genes that do not meet the threshold from the data set
#expr <- expr[which(expressed),]
#meta_genes <- meta_genes[which(expressed),]

# (Alternative) Adding a column at the gene metadata information table to mark whether it is expressed
meta_genes$expressed <- expressed

# ====================================================================================================
# ### Correlation clustering
# - **Pearson** correlation on log TPM (linear) and **Spearman** correlation (rank-based, robust to skewed expression).
# - Hierarchical clustering with 1 − correlation as distance.
# ====================================================================================================

# 2 · Correlation clustering of all samples: Pearson on log(TPM + 1), Spearman on TPM;
# distance = 1 - correlation. The Spearman tree is then relabelled by group and by dataset.
options(repr.plot.width = 16, repr.plot.height = 7)

corr_pearson <- cor(log1p(expr[meta_genes$expressed,]))
corr_spearman <- cor(expr[meta_genes$expressed,], method = "spearman")

# Apply hierarchical clustering to the samples, given the correlation coefficients as similarities
hcl_pearson <- hclust(as.dist(1 - corr_pearson))
hcl_spearman <- hclust(as.dist(1 - corr_spearman))

layout(matrix(1:2,nrow=1))
plot(hcl_pearson)
plot(hcl_spearman)



# Plotting the meta data labels instead of the accessions
layout(matrix(1:2,nrow=1))
plot(hcl_spearman, labels = meta$group)
plot(hcl_spearman, labels = meta$dataset)

# **Result:** samples cluster by dataset, a strong batch effect between the two BioProjects. Datasets are therefore clustered separately and compared only through fold changes.

# 2 · The same clustering within each dataset (the batch effect hides group structure otherwise);
# labels = group + sex
options(repr.plot.width = 16, repr.plot.height = 7)

for (ds in unique(meta$dataset)) {
  # Samples of this dataset
  s <- meta$SRR[meta$dataset == ds]
  x <- expr[meta_genes$expressed, s]
  x <- x[apply(x, 1, var) > 0, ]   # dropping genes that are constant in this dataset

  # Calculating PCC and SCC between samples
  corr_pearson_ds  <- cor(log1p(x))
  corr_spearman_ds <- cor(x, method = "spearman")

  # Hierarchical clustering
  hcl_pearson_ds  <- hclust(as.dist(1 - corr_pearson_ds))
  hcl_spearman_ds <- hclust(as.dist(1 - corr_spearman_ds))

  # Labels: group and sex of each sample
  labs <- paste(meta$group, meta$sex)[match(s, meta$SRR)]

  layout(matrix(1:2, nrow = 1))
  par(mar = c(10, 4, 3, 1))
  plot(hcl_pearson_ds,  labels = labs, hang = -1, main = paste("Pearson -", ds),  xlab = "", sub = "")
  plot(hcl_spearman_ds, labels = labs, hang = -1, main = paste("Spearman -", ds), xlab = "", sub = "")
}

# **Result:** within each dataset, neither clustering consistently separates age or intervention groups.

# ====================================================================================================
# ### PCA
# Log TPM, scaled; all samples and each dataset separately.
# ====================================================================================================

# 2 · PCA on log(TPM + 1) of expressed genes, scaled (every gene weighted equally);
# all samples, then each dataset. Function: PCA plot for a subset of samples
plot_pca <- function(samples, title) {
  x <- log1p(t(expr[meta_genes$expressed, samples]))
  x <- x[, apply(x, 2, var) > 0]      # remove genes that are constant in this subset
  pca <- prcomp(x, center = TRUE, scale. = TRUE)
  var_expl <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  ggplot(data.frame(pca$x, meta[match(samples, meta$SRR), ])) +
    geom_point(aes(x = PC1, y = PC2, color = group, shape = cohort), size = 4) +
    labs(title = title,
         x = paste0("PC1 (", var_expl[1], "%)"),
         y = paste0("PC2 (", var_expl[2], "%)")) +
    theme_minimal()
}

plot_pca(meta$SRR, "All samples")
plot_pca(meta$SRR[meta$dataset == "PRJNA875066"], "Aging dataset")
plot_pca(meta$SRR[meta$dataset == "PRJNA946653"], "Intervention dataset")

# **Result:** across all samples, PC1 separates the two cohorts (the batch effect). Within each dataset, PC1 separates neither age groups nor interventions: the biological effects are subtle relative to differences between mice.

# ====================================================================================================
# ## 3 · Differential expression (DESeq2)
# - Input: RSEM expected counts, rounded to integers.
# - Genes kept: ≥ 10 reads in ≥ 3 samples.
# - DEG: Benjamini–Hochberg padj < 0.05, used throughout.
#
# | Comparison | Groups |
# |---|---|
# | Aging | 21 vs 18 months (to justify pooling), then old (18 + 21) vs young (3 months) |
# | Dietary restriction | aDR vs ad libitum (AL) |
# | Young plasma | YMP vs PBS |
# ====================================================================================================

# 3 · Differential expression (DESeq2, Wald test), genes with >= 10 counts in >= 3 samples. Builds:
#   de_age18, de_age21, de_age21v18 (three age groups; only to decide on pooling)
#   de_old (old 18 + 21 mo vs young 3 mo), de_diet (aDR vs AL), de_plasma (YMP vs PBS)
library(DESeq2)

# Function: extract one contrast and add gene symbols
get_res <- function(dds, contrast) {
  # alpha = 0.05 matches the padj < 0.05 cutoff used everywhere (DESeq2's default is 0.1)
  res <- as.data.frame(results(dds, contrast = contrast, alpha = 0.05))
  res$gene_id <- rownames(res)
  res$symbol  <- meta_genes$symbol[match(res$gene_id, meta_genes$gene_id)]
  res[order(res$padj), ]
}

# Function: fit DESeq2 for one intervention, "treat" vs "ctrl"
# (each intervention cohort is single-sex: diet = female, plasma = male, so no sex term is needed)
# The fitted objects are reused later (shrinkage, camera), so each model is fitted only once
fit_interv <- function(samples, treat_group) {
  cd <- meta[match(samples, meta$SRR), ]
  cd$condition <- factor(ifelse(cd$group == treat_group, "treat", "ctrl"), levels = c("ctrl", "treat"))
  cts <- round(counts[, samples])        # RSEM expected counts are not integers
  cts <- cts[rowSums(cts >= 10) >= 3, ]  # keep genes with at least 10 reads in at least 3 samples
  dds <- DESeqDataSetFromMatrix(cts, colData = cd, design = ~ condition)
  # DESeq() does three things:
  #   1. normalisation: size factors correct for different sequencing depths
  #   2. dispersion: how much each gene varies between replicates; with 3-5 mice per group, information is shared across genes
  #   3. test: a model per gene, Wald test of the condition effect
  DESeq(dds)
}

# Aging samples, genes with at least 10 reads in at least 3 samples
s_age <- meta$SRR[meta$cohort == "Aging"]
cd_age <- meta[match(s_age, meta$SRR), ]
cts_age <- round(counts[, s_age])
cts_age <- cts_age[rowSums(cts_age >= 10) >= 3, ]

# Step 1 - aging with the three age groups kept separate
# (used only to test whether 18 and 21 months differ)
cd_age$age <- factor(paste0("m", cd_age$age_months), levels = c("m3", "m18", "m21"))
dds_age <- cached("dds_age",
                  DESeq(DESeqDataSetFromMatrix(cts_age, colData = cd_age, design = ~ age)))

de_age18    <- get_res(dds_age, c("age", "m18", "m3"))   # 18 mo vs 3 mo
de_age21    <- get_res(dds_age, c("age", "m21", "m3"))   # 21 mo vs 3 mo
de_age21v18 <- get_res(dds_age, c("age", "m21", "m18"))  # 21 mo vs 18 mo

# Step 2 - 18 and 21 months binned into one "old" group vs 3 months ("young")
# (used for the rest of the analysis)
cd_age$age_bin <- factor(ifelse(cd_age$age_months == 3, "young", "old"), levels = c("young", "old"))
dds_old <- cached("dds_old",
                  DESeq(DESeqDataSetFromMatrix(cts_age, colData = cd_age, design = ~ age_bin)))

de_old <- get_res(dds_old, c("age_bin", "old", "young"))   # old (18 + 21 mo) vs young (3 mo)

# Dietary restriction: aDR vs AL
s_diet   <- meta$SRR[meta$group %in% c("AL", "aDR")]
dds_diet <- cached("dds_diet", fit_interv(s_diet, "aDR"))
de_diet  <- get_res(dds_diet, c("condition", "treat", "ctrl"))

# Young plasma: YMP vs PBS
s_plasma   <- meta$SRR[meta$group %in% c("PBS", "YMP")]
dds_plasma <- cached("dds_plasma", fit_interv(s_plasma, "YMP"))
de_plasma  <- get_res(dds_plasma, c("condition", "treat", "ctrl"))

# Number of significant genes per comparison
sapply(list(age18 = de_age18, age21 = de_age21, age21v18 = de_age21v18, old = de_old,
            diet = de_diet, plasma = de_plasma),
       function(r) sum(r$padj < 0.05, na.rm = TRUE))

head(de_old[, c("symbol", "log2FoldChange", "padj")], 20)

# ====================================================================================================
# ### Pooling the 18- and 21-month groups
# The two old groups are pooled if:
# 1. few genes differ between 21 and 18 months, compared with each old group vs 3 months, and
# 2. their fold changes vs 3 months are strongly correlated and close to the identity line.
#
# All later analyses use the pooled comparison, old vs young (`de_old`).
# ====================================================================================================

# 3 · Pool 18 and 21 months? (a) DEGs 21 vs 18 compared with each vs 3 months;
# (b) Spearman correlation and scatter of their fold changes vs 3 months (dashed line = identical change)
options(repr.plot.width = 7, repr.plot.height = 6)

# (a) Number of DEGs: 21 vs 18 mo next to each old group vs 3 mo
n_deg <- function(r) c(padj_0.05            = sum(r$padj < 0.05, na.rm = TRUE),
                       padj_0.05_and_2_fold = sum(r$padj < 0.05 & abs(r$log2FoldChange) > 1, na.rm = TRUE))
sapply(list("18 vs 3 mo" = de_age18, "21 vs 3 mo" = de_age21, "21 vs 18 mo" = de_age21v18), n_deg)

# (b) Do 18 and 21 mo show the same aging changes? Fold changes vs 3 mo
both <- inner_join(de_age18 %>% select(gene_id, symbol, lfc18 = log2FoldChange, padj18 = padj),
                   de_age21 %>% select(gene_id, lfc21 = log2FoldChange, padj21 = padj),
                   by = "gene_id")
sig <- both %>% filter(padj18 < 0.05 | padj21 < 0.05)   # genes that change with age in either group

c(rho_all_genes   = cor(both$lfc18, both$lfc21, method = "spearman", use = "complete.obs"),
  rho_aging_genes = cor(sig$lfc18,  sig$lfc21,  method = "spearman", use = "complete.obs"))

ggplot(sig, aes(x = lfc18, y = lfc21)) +
  geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.6) +
  labs(title = "Aging changes at 18 and 21 months",
       subtitle = "Genes significant (padj < 0.05) in 18 vs 3 or 21 vs 3 mo; dashed line = identical change",
       x = "log2FC 18 vs 3 mo", y = "log2FC 21 vs 3 mo") +
  theme_minimal()

# 3 · Volcano plots: coloured = padj < 0.05 and |log2FC| > 1; the 15 most significant are labelled
library(ggrepel)
options(repr.plot.width = 10, repr.plot.height = 7)

# Function: volcano plot of one DE result table
plot_volcano <- function(res, title, padj_cut = 0.05, lfc_cut = 1, n_labels = 15) {
  res <- res[!is.na(res$padj), ]
  res$direction <- "n.s."
  res$direction[res$padj < padj_cut & res$log2FoldChange >  lfc_cut] <- "up"
  res$direction[res$padj < padj_cut & res$log2FoldChange < -lfc_cut] <- "down"
  top <- head(res[res$direction != "n.s.", ], n_labels)   # res is sorted by padj

  ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = direction)) +
    geom_point(size = 1, alpha = 0.6) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", color = "grey50") +
    geom_text_repel(data = top, aes(label = symbol), size = 3, max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(values = c(down = "#2166ac", n.s. = "grey75", up = "#b2182b")) +
    labs(title = title,
         subtitle = paste0(sum(res$direction == "up"), " up, ",
                           sum(res$direction == "down"), " down"),
         x = "log2 fold change", y = "-log10 adjusted p-value") +
    theme_minimal()
}

plot_volcano(de_age21v18, "Aging: 21 mo vs 18 mo (step 1)")
plot_volcano(de_old,      "Aging: old (18 + 21 mo) vs young (3 mo)")
plot_volcano(de_diet,     "Dietary restriction: aDR vs AL")
plot_volcano(de_plasma,   "Young plasma: YMP vs PBS")

# **Volcano plots:** dashed lines at padj = 0.05 and |log2FC| = 1 (two-fold). Genes passing both are coloured (red up, blue down); the 15 most significant are labelled.

# 3 · Heatmaps of the top 50 DEGs per comparison: log(TPM + 1), z-scored per gene; bar = sample group
library(gplots)
options(repr.plot.width = 10, repr.plot.height = 12)

# One colour per experimental group (for the bar above the heatmap)
group_cols <- c("aging 3mo" = "#c7e9c0", "aging 18mo" = "#74c476", "aging 21mo" = "#238b45",
                "AL" = "#fdd0a2", "aDR" = "#e6550d",
                "PBS" = "#c6dbef", "YMP" = "#3182bd")

# Function: heatmap of the top n significant genes of one DE result, across the chosen samples
plot_heatmap <- function(res, samples, title, n = 50, padj_cut = 0.05) {
  top <- head(res[!is.na(res$padj) & res$padj < padj_cut, ], n)
  if (nrow(top) < 2) {
    message(title, ": fewer than 2 significant genes, no heatmap")
    return(invisible(NULL))
  }

  # log TPM, then z-score per gene (row): how much higher/lower than the gene's mean
  x <- log1p(expr[top$gene_id, samples])
  z <- t(scale(t(x)))
  rownames(z) <- ifelse(is.na(top$symbol), top$gene_id, top$symbol)

  grp <- as.character(meta$group[match(samples, meta$SRR)])

  heatmap.2(z,
            col = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
            breaks = seq(-2.5, 2.5, length.out = 101),   # cap colours at |z| = 2.5
            scale = "none", trace = "none", density.info = "none",
            ColSideColors = group_cols[grp],
            labCol = grp, cexCol = 0.9, cexRow = 0.7,
            margins = c(8, 8), key.title = NA, key.xlab = "z-score",
            main = title)
  legend("topright", legend = unique(grp), fill = group_cols[unique(grp)],
         bty = "n", cex = 0.8, inset = c(0, 0))
}

s_aging <- meta$SRR[meta$cohort == "Aging"]

plot_heatmap(de_old,    s_aging,  "Top aging genes (old vs young), all aging samples")
plot_heatmap(de_diet,   s_diet,   "Top genes aDR vs AL")
plot_heatmap(de_plasma, s_plasma, "Top genes YMP vs PBS")

# **Result: aging signature.** Mostly upregulated genes, forming a coherent signature:
#
# | Category | Genes |
# |---|---|
# | Microglial activation (disease-associated microglia) | *Clec7a, Cst7, Trem2, Tyrobp, Lpl, Gpnmb, Lgals3, Cd68, Lilrb4a, Fcgr2b* |
# | Complement system | *C4b, C1qa, C1qb, C1qc* |
# | MHC class I / antigen presentation | *H2-Q6, H2-Q7* |
# | Lysosomal proteases | *Ctss, Ctsd, Ctsz* |
# | Reactive astrocytes | *Gfap* |
#
# - Together: the neuroinflammatory ("inflammaging") signature of the aging brain.
# - Down with age: oligodendrocyte / myelin genes (*Mal, Cldn11, Gjc2*) and proliferation genes (*Mcm3, Ccnd1*), in line with the white-matter vulnerability reported by Hahn et al. (2023).
# - Few genes differ between 21 and 18 months, and the two groups intermix in the clustering: most changes are established by 18 months, which supports pooling.

# **Summary of differential expression**
# - **Aging:** neuroinflammation (microglial activation, complement, MHC class I) and declining myelin genes; largely established by 18 months.
# - **Dietary restriction:** robust but modest; dominated by circadian, glucocorticoid-responsive and metabolic genes, not by reversal of inflammatory genes.
# - **Young plasma:** few significant genes; may partly reflect differences in tissue composition.
#
# None of the top aging genes (e.g. *C4b, Clec7a, Trem2, C1q*) are among the top intervention genes, so neither intervention strongly reverses the core aging signature at the single-gene level. Section 5 tests this across all aging genes.

# ====================================================================================================
# ## 4 · Gene set enrichment (GSEA)
# GSEA uses all genes, ranked, not only significant ones, so it detects coordinated modest shifts across many genes (as in the dietary restriction response).
#
# - **Gene sets:** mouse-native MSigDB, Hallmark (MH, 50 sets) and GO Biological Process (M5, GO:BP).
# ====================================================================================================

# - **Ranking:** DESeq2 Wald statistic (combines effect size and its uncertainty).
# - **Test:** fgsea, sets of 15–500 genes; are a set's genes concentrated at the top (up) or bottom (down) of the ranking?
# - **NES** (normalized enrichment score): positive = up in old or treated mice, negative = down.
# - **padj:** Benjamini–Hochberg across gene sets; significant at padj < 0.05.

# 4 · GSEA: mouse MSigDB gene sets (Hallmark, GO BP); all genes ranked by the DESeq2 Wald statistic;
# fgsea for each comparison (old vs young, aDR, YMP). Results are cached.
library(fgsea)
library(msigdbr)

# Gene sets for mouse: Hallmark (50 broad processes) and GO Biological Process
# db_species = "MM" uses the mouse-native MSigDB (MH, M5) instead of human sets mapped to mouse genes
hallmark_sets <- cached("hallmark_sets", {
  hallmark <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH")
  split(hallmark$gene_symbol, hallmark$gs_name)
})
gobp_sets <- cached("gobp_sets", {
  gobp <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "M5", subcollection = "GO:BP")
  split(gobp$gene_symbol, gobp$gs_name)
})

# Function: ranked gene list from a DESeq2 result (Wald statistic, one value per gene symbol)
make_ranks <- function(res) {
  res <- res[!is.na(res$stat) & !is.na(res$symbol), ]
  res <- res[order(-abs(res$stat)), ]
  res <- res[!duplicated(res$symbol), ]           # keep the strongest entry if a symbol appears twice
  sort(setNames(res$stat, res$symbol), decreasing = TRUE)
}

# Function: run GSEA for one comparison
run_gsea <- function(res, sets) {
  set.seed(42)                                    # GSEA uses permutations, so fix the random seed
  out <- fgsea(pathways = sets, stats = make_ranks(res), minSize = 15, maxSize = 500)
  out[order(out$padj), ]
}

de_list <- list(old = de_old, diet = de_diet, plasma = de_plasma)

gsea_hallmark <- cached("gsea_hallmark", lapply(de_list, run_gsea, sets = hallmark_sets))
gsea_gobp     <- cached("gsea_gobp",     lapply(de_list, run_gsea, sets = gobp_sets))

# Number of significant gene sets per comparison
sapply(gsea_hallmark, function(r) sum(r$padj < 0.05))
sapply(gsea_gobp,     function(r) sum(r$padj < 0.05))

# 4 · GSEA figures: (1) Hallmark dot plot across the three comparisons (colour = NES, size = -log10 padj);
# (2) the 10 most up- and downregulated significant GO BP terms per comparison
options(repr.plot.width = 12, repr.plot.height = 10)

# 1) Hallmark overview: NES of all comparisons side by side
nes <- do.call(rbind, lapply(names(gsea_hallmark), function(n) {
  data.frame(comparison = n, pathway = gsea_hallmark[[n]]$pathway,
             NES = gsea_hallmark[[n]]$NES, padj = gsea_hallmark[[n]]$padj)
}))
nes$pathway    <- sub("HALLMARK_", "", nes$pathway)
nes$comparison <- factor(nes$comparison, levels = names(de_list))

# keep pathways significant in at least one comparison
keep <- unique(nes$pathway[nes$padj < 0.05])

ggplot(nes[nes$pathway %in% keep, ], aes(x = comparison, y = reorder(pathway, NES))) +
  geom_point(aes(size = -log10(padj), color = NES)) +
  scale_color_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
  labs(title = "Hallmark gene sets (GSEA)", x = NULL, y = NULL,
       size = "-log10 padj", color = "NES") +
  theme_minimal()

# 2) Top GO terms per comparison
plot_top_go <- function(gs, title, n = 10) {
  gs <- gs[gs$padj < 0.05, ]
  top <- rbind(head(gs[order(-gs$NES), ], n), head(gs[order(gs$NES), ], n))
  top <- top[!duplicated(top$pathway) & top$padj < 0.05, ]
  if (nrow(top) == 0) { message(title, ": no significant GO terms"); return(invisible(NULL)) }
  top$pathway <- tolower(gsub("_", " ", sub("GOBP_", "", top$pathway)))
  ggplot(top, aes(x = NES, y = reorder(pathway, NES), fill = NES > 0)) +
    geom_col() +
    scale_fill_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "#2166ac"),
                      labels = c(`TRUE` = "up", `FALSE` = "down"), name = NULL) +
    labs(title = title, x = "Normalized enrichment score (NES)", y = NULL) +
    theme_minimal()
}

plot_top_go(gsea_gobp$old,    "GO BP: old (18 + 21 mo) vs young (3 mo)")
plot_top_go(gsea_gobp$diet,   "GO BP: aDR vs AL")
plot_top_go(gsea_gobp$plasma, "GO BP: YMP vs PBS")

# **Figures:** Hallmark overview, colour = NES and size = −log10 padj, only sets significant in at least one comparison. GO bar plots: the 10 most up- and downregulated significant terms per comparison.

# ====================================================================================================
# ## 5 · Intervention effects vs aging
# The datasets come from separate experiments, so only fold changes are compared, never expression values.
#
# ### 5a · Version 1 (superseded by 5b, kept for reference)
# Known issues:
# 1. Raw (unshrunken) log2 fold changes, noisy for low-count genes.
# 2. "Reversed" based on the sign of the fold change alone, without evidence from the intervention.
# 3. Binomial and correlation tests treat genes as independent, so p-values are too small.
# 4. GO enrichment of reversed genes uses all genes as background, so it mostly reproduces the aging signature.
# ====================================================================================================

# 5a (Version 1, superseded by 5b): aging vs intervention log2FC (unshrunken) over the aging DEGs:
# Spearman correlation, % with opposite sign (binomial test), and the same correlation over all genes
# Join log2 fold changes of all comparisons per gene
lfc <- de_old %>% select(gene_id, symbol, old = log2FoldChange, padj_old = padj) %>%
  inner_join(de_diet    %>% select(gene_id, diet   = log2FoldChange), by = "gene_id") %>%
  inner_join(de_plasma  %>% select(gene_id, plasma = log2FoldChange), by = "gene_id")

# Aging genes: significant in old (18 + 21 mo) vs young (3 mo)
age_genes <- lfc %>% filter(padj_old < 0.05)
nrow(age_genes)

# Function: summary of how an intervention relates to an aging signature
compare_to_aging <- function(age, interv) {
  ct <- cor.test(age, interv, method = "spearman", exact = FALSE)
  reversed <- sum(sign(age) != sign(interv))
  bt <- binom.test(reversed, length(age), p = 0.5)       # 50 % expected by chance
  data.frame(rho = round(ct$estimate, 3), p_cor = signif(ct$p.value, 3),
             n_genes = length(age), n_reversed = reversed,
             pct_reversed = round(100 * reversed / length(age), 1),
             p_binom = signif(bt$p.value, 3))
}

summary_tab <- rbind(
  "diet vs aging"   = compare_to_aging(age_genes$old, age_genes$diet),
  "plasma vs aging" = compare_to_aging(age_genes$old, age_genes$plasma)
)
summary_tab

# Same correlation on all genes, as a background check
c(diet_all   = cor(lfc$old, lfc$diet,   method = "spearman"),
  plasma_all = cor(lfc$old, lfc$plasma, method = "spearman"))

# **Metrics** (aging genes only: old vs young, padj < 0.05):
# - **rho:** Spearman correlation of aging and intervention log2FC. Negative = rejuvenating, positive = mimics aging, ≈ 0 = no relationship.
# - **pct_reversed / p_binom:** % of aging genes whose direction the intervention reverses; two-sided binomial test against 50 %.
# - **All-genes correlation:** background control; if as strong as over aging genes, the effect is not specific to aging.

# 5a (superseded): scatter of aging vs intervention log2FC per aging gene; green = opposite sign, orange = same
options(repr.plot.width = 12, repr.plot.height = 6)

# Long format: one row per gene x intervention
plot_df <- rbind(
  data.frame(age_genes[, c("gene_id", "symbol")], interv = "aDR vs AL",
             age_lfc = age_genes$old, int_lfc = age_genes$diet),
  data.frame(age_genes[, c("gene_id", "symbol")], interv = "YMP vs PBS",
             age_lfc = age_genes$old, int_lfc = age_genes$plasma)
)
plot_df$effect <- ifelse(sign(plot_df$age_lfc) != sign(plot_df$int_lfc), "reversed", "same direction")

# Label the genes with the strongest reversal in each panel
lab_df <- plot_df %>% filter(effect == "reversed") %>%
  group_by(interv) %>% slice_max(abs(age_lfc * int_lfc), n = 8)

ggplot(plot_df, aes(x = age_lfc, y = int_lfc)) +
  geom_hline(yintercept = 0, color = "grey50") + geom_vline(xintercept = 0, color = "grey50") +
  geom_point(aes(color = effect), alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
  geom_text_repel(data = lab_df, aes(label = symbol), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c(reversed = "#1b9e77", `same direction` = "#d95f02")) +
  facet_wrap(~ interv) +
  labs(x = "log2FC aging (old vs young)", y = "log2FC intervention", color = NULL,
       title = "Do interventions reverse aging-associated gene expression?") +
  theme_minimal()

# **Figure:** each point is an aging gene; x = log2FC with age, y = log2FC with the intervention. Green (upper-left, lower-right) = reversed; orange = same direction as aging. Black line: linear fit with 95 % CI (negative slope = rejuvenating trend). The 8 most strongly reversed genes per panel are labelled.

# 5a (superseded): the same summary for genes up and for genes down with age; lists of sign-reversed genes
# Split aging genes by direction of the aging change (old vs young)
age_up   <- age_genes %>% filter(old > 0)
age_down <- age_genes %>% filter(old < 0)

split_tab <- rbind(
  "diet, genes UP with age"     = compare_to_aging(age_up$old,   age_up$diet),
  "diet, genes DOWN with age"   = compare_to_aging(age_down$old, age_down$diet),
  "plasma, genes UP with age"   = compare_to_aging(age_up$old,   age_up$plasma),
  "plasma, genes DOWN with age" = compare_to_aging(age_down$old, age_down$plasma)
)
split_tab

# Reversed genes in each half: the intervention changes them in the opposite direction to aging
get_reversed <- function(df, interv) {
  g <- df$symbol[sign(df$old) != sign(df[[interv]])]
  unique(na.omit(g))
}

reversed_genes <- list(
  "diet, genes UP with age"     = get_reversed(age_up,   "diet"),
  "diet, genes DOWN with age"   = get_reversed(age_down, "diet"),
  "plasma, genes UP with age"   = get_reversed(age_up,   "plasma"),
  "plasma, genes DOWN with age" = get_reversed(age_down, "plasma")
)
sapply(reversed_genes, length)

# 5a (superseded): save the split table as CSV and PNG
# Save as table (opens in Excel)
write.csv(signif(split_tab, 3), "split_table.csv")

# Save as picture
library(gridExtra)
png("split_table.png", width = 1200, height = 250, res = 120)
grid.table(signif(split_tab, 3))
dev.off()

# **Result** (written for the first run; the re-run with pooled aging groups no longer supports all of it):
# - Young plasma reverses both parts of the aging signature: it attenuates age-induced (largely inflammatory) genes and, more consistently, restores age-repressed genes.
# - Dietary restriction partly attenuates age-induced inflammatory genes but reinforces the decline of age-repressed genes, mostly myelin and oligodendrocyte genes.

# 5a (superseded): GO over-representation of the sign-reversed genes
# Issue: the background (universe) is all genes, but reversed genes can only come from the aging genes
universe <- unique(na.omit(lfc$symbol))   # all genes tested in all comparisons

ora_reversed <- lapply(reversed_genes, function(g)
  fora(pathways = gobp_sets, genes = g, universe = universe, minSize = 10))

# Top 10 GO terms per list
lapply(ora_reversed, function(r) head(r[order(r$padj), c("pathway", "overlap", "size", "padj")], 10))

# 5a (superseded): dot plot of the top 10 GO terms per list, saved as reversed_genes_GO.png
options(repr.plot.width = 14, repr.plot.height = 10)

# Combine the top 10 terms of each list into one table
ora_df <- do.call(rbind, lapply(names(ora_reversed), function(n) {
  r <- ora_reversed[[n]]
  r <- head(r[order(r$padj), ], 10)
  data.frame(list = n,
             term = tolower(gsub("_", " ", sub("GOBP_", "", r$pathway))),
             overlap = r$overlap,
             ratio = r$overlap / length(reversed_genes[[n]]),   # share of the list in this term
             padj = r$padj)
}))

ggplot(ora_df, aes(x = ratio, y = reorder(term, ratio))) +
  geom_point(aes(size = overlap, color = -log10(padj))) +
  scale_color_viridis_c() +
  facet_wrap(~ list, scales = "free_y", ncol = 2) +
  labs(title = "GO terms enriched among reversed aging genes",
       x = "Gene ratio (genes in term / genes in list)", y = NULL,
       size = "Genes in term", color = "-log10 padj") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

# Save the plot
ggsave("reversed_genes_GO.png", width = 14, height = 10, dpi = 150)

# ====================================================================================================
# ### 5b · Version 2
# **Question:** does each intervention shift the cerebellar transcriptome back towards a young state?
#
# **Design:** the aging cohort has both sexes at 3, 18 and 21 months. The intervention cohorts are 19-month-old mice of one sex (aDR female, YMP male). Their controls (AL, PBS) are old too, so intervention fold changes describe changes within old animals and can be compared with aging fold changes.
#
# **Levels, from most to least robust:**
# 1. **Preparation:** shrunken fold changes on a common gene universe; aging signature split by direction.
# 2. **Signature (main test):** are the aging-up and aging-down genes, as sets, reversed?
# 3. **Pathway:** which Hallmark and GO processes are reversed or reinforced?
# 4. **Gene:** which aging genes are reversed, and what are their functions?
#
# Three controls: specificity for aging genes, robustness to the aging-signature definition, cell-type composition.
#
# | Issue in Version 1 | Change in Version 2 |
# |---|---|
# | Raw, noisy fold changes | apeglm-shrunken fold changes |
# | Reversal from the sign alone | camera/GSEA on gene sets; at gene level, intervention padj < 0.05 required |
# | Genes assumed independent | camera accounts for inter-gene correlation; specificity tested against random gene sets |
# | GO background = all genes | GO background = aging genes of the same direction |
# | Sex not considered | Comparison with same-sex aging signatures |
# ====================================================================================================

# ====================================================================================================
# #### Step 1 · Shrunken fold changes and aging signature
# With 3–5 mice per group, raw log2FC of low-count or variable genes are noisy and can be extreme.
# - **Models:** the aDR vs AL and YMP vs PBS models from section 3 (single-sex cohorts, so no sex term).
# - **Shrinkage:** apeglm (`lfcShrink(type = "apeglm")`; Zhu, Ibrahim & Love, 2019) pulls uncertain fold changes towards zero; well-measured genes keep theirs. p-values, and so the significant genes, are unchanged.
# - **Common gene universe:** genes passing the filter in all three comparisons.
# - **Aging signature:** old vs young padj < 0.05, split into **aging_up** (higher in old) and **aging_down** (lower in old).
# ====================================================================================================

# 5b step 1: apeglm-shrunken log2FC for the three comparisons, on genes tested in all three (common_genes);
# aging signature = old vs young padj < 0.05, split into aging_up / aging_down
library(apeglm)

# dds_old, dds_diet and dds_plasma are the models already fitted in the DESeq2 cell above

# Genes that pass the expression filter in all three datasets
common_genes <- Reduce(intersect, list(rownames(dds_old), rownames(dds_diet), rownames(dds_plasma)))

# Function: apeglm-shrunken log2 fold changes, named by gene_id
shrink <- function(dds, coef) {
  r <- as.data.frame(lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE))
  setNames(r$log2FoldChange, rownames(r))
}
shr_old    <- cached("shr_old",    shrink(dds_old,    "age_bin_old_vs_young"))
shr_diet   <- cached("shr_diet",   shrink(dds_diet,   "condition_treat_vs_ctrl"))
shr_plasma <- cached("shr_plasma", shrink(dds_plasma, "condition_treat_vs_ctrl"))

# One table: shrunken fold changes plus the (unchanged) p-values of each comparison
v2 <- data.frame(gene_id  = common_genes,
                 symbol   = meta_genes$symbol[match(common_genes, meta_genes$gene_id)],
                 old      = shr_old[common_genes],
                 diet     = shr_diet[common_genes],
                 plasma   = shr_plasma[common_genes],
                 padj_old = de_old$padj[match(common_genes, de_old$gene_id)],
                 padj_diet   = de_diet$padj[match(common_genes, de_diet$gene_id)],
                 padj_plasma = de_plasma$padj[match(common_genes, de_plasma$gene_id)],
                 baseMean = de_old$baseMean[match(common_genes, de_old$gene_id)])

# Aging signature: old vs young padj < 0.05, split by direction
aging_sets <- list(aging_up   = v2$gene_id[which(v2$padj_old < 0.05 & v2$old > 0)],
                   aging_down = v2$gene_id[which(v2$padj_old < 0.05 & v2$old < 0)])

c(common_genes = length(common_genes), lengths(aging_sets))

# ====================================================================================================
# #### Step 2 · Signature level: is the aging signature reversed as a whole?
# The main test: aging_up and aging_down are treated as two gene sets and tested in the intervention data, giving one answer per intervention and direction.
# - **camera** (limma; Wu & Smyth, 2012) on voom-transformed intervention counts: are the set's genes, on average, more up- or downregulated than other genes? Accounts for inter-gene correlation. padj = Benjamini–Hochberg across the two sets.
# - **GSEA** (fgsea) of the two sets on the intervention ranking (Wald statistic), as a complementary rank-based test with enrichment plots.
#
# | Gene set | Rejuvenating | Mimics aging |
# |---|---|---|
# | aging_up | "Down" / negative NES | "Up" / positive NES |
# | aging_down | "Up" / positive NES | "Down" / negative NES |
#
# **Enrichment plots:** genes ordered from most up- (left) to most downregulated (right) by the intervention; ticks = aging genes; green curve = running enrichment score. A peak above zero = set concentrated among upregulated genes; a dip below zero = among downregulated genes.
# ====================================================================================================

# 5b step 2 (main test): is the aging signature reversed as a whole?
# camera (on voom-transformed counts) and fgsea test aging_up / aging_down in each intervention;
# then enrichment plots and a bubble plot
library(limma)
library(edgeR)
library(gridExtra)
options(repr.plot.width = 12, repr.plot.height = 8)

# Function: camera test of gene sets in one intervention dataset
# camera accounts for correlation between genes, so its p-values are not inflated like a per-gene test
run_camera <- function(dds, sets) {
  sets <- sets[lengths(sets) >= 5]                       # camera needs a few genes per set
  y <- calcNormFactors(DGEList(counts(dds)[common_genes, ]))
  design <- model.matrix(~ condition, data = as.data.frame(colData(dds)))
  v <- voom(y, design)
  res <- camera(v, ids2indices(sets, rownames(v)), design, contrast = "conditiontreat")
  res$set <- rownames(res)
  res$padj <- p.adjust(res$PValue, method = "BH")      # adjusted across the tested sets
  # Rejuvenating = aging-up genes go down, or aging-down genes go up
  res$rejuvenating <- (res$set == "aging_up" & res$Direction == "Down") |
                      (res$set == "aging_down" & res$Direction == "Up")
  res[, c("set", "NGenes", "Direction", "padj", "rejuvenating")]
}

camera_main <- rbind(cbind(intervention = "aDR vs AL",  run_camera(dds_diet,   aging_sets)),
                     cbind(intervention = "YMP vs PBS", run_camera(dds_plasma, aging_sets)))
rownames(camera_main) <- NULL
camera_main

# Same question with GSEA: where do the aging genes sit in the intervention ranking?
ranks_id <- function(res) {
  r <- res[res$gene_id %in% common_genes & !is.na(res$stat), ]
  sort(setNames(r$stat, r$gene_id), decreasing = TRUE)
}
set.seed(42)
gsea_sig <- lapply(list(diet = de_diet, plasma = de_plasma), function(r) fgsea(aging_sets, ranks_id(r)))
lapply(gsea_sig, function(g) g[, c("pathway", "size", "NES", "padj")])

# Enrichment plots: rejuvenation = aging_up curve below zero, aging_down curve above zero
plots <- list()
for (i in names(gsea_sig)) for (s in names(aging_sets)) {
  g <- gsea_sig[[i]][gsea_sig[[i]]$pathway == s, ]
  plots[[paste(i, s)]] <- plotEnrichment(aging_sets[[s]], ranks_id(get(paste0("de_", i)))) +
    labs(title = paste0(i, ": ", s), subtitle = sprintf("NES = %.2f, padj = %.2g", g$NES, g$padj))
}
grid.arrange(grobs = plots, ncol = 2)

# Bubble plot: one dot per intervention x aging set
# fill = NES (blue = lower after the intervention, red = higher), size = -log10 padj, black outline = padj < 0.05
# Rejuvenation = aging_up genes blue (pushed down), aging_down genes red (pushed up)
options(repr.plot.width = 8, repr.plot.height = 5)

bubble <- do.call(rbind, lapply(names(gsea_sig), function(i) {
  g <- as.data.frame(gsea_sig[[i]])
  data.frame(intervention = c(diet = "aDR vs AL", plasma = "YMP vs PBS")[[i]],
             set = g$pathway, n_genes = g$size, NES = g$NES, padj = g$padj)
}))
bubble$effect <- ifelse((bubble$set == "aging_up" & bubble$NES < 0) |
                        (bubble$set == "aging_down" & bubble$NES > 0), "rejuvenating", "aging-like")
bubble$set <- factor(bubble$set, levels = c("aging_down", "aging_up"),
                     labels = c("aging_down\n(lower in old mice)", "aging_up\n(higher in old mice)"))
lim <- max(abs(bubble$NES))   # symmetric colour scale around 0

ggplot(bubble, aes(x = intervention, y = set)) +
  geom_point(aes(size = -log10(padj), fill = NES, color = padj < 0.05), shape = 21, stroke = 1.2) +
  geom_text(aes(label = sprintf("NES = %.2f, padj = %.2g\n%s", NES, padj, effect)),
            nudge_y = -0.35, size = 3.2) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-lim, lim)) +
  scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "grey70"), name = "padj < 0.05") +
  scale_size_continuous(range = c(6, 18), name = "-log10 padj") +
  scale_y_discrete(expand = expansion(add = 0.7)) +
  labs(title = "Is the aging signature reversed by the interventions? (GSEA)",
       subtitle = "Rejuvenation: aging_up blue (pushed down), aging_down red (pushed up)",
       x = NULL, y = NULL) +
  theme_minimal()

# ====================================================================================================
# #### Step 3 · Pathway level: which processes are reversed?
# Pathway scores aggregate many genes, so they are less noisy than single genes. For each Hallmark and GO BP set, the aging NES is compared with the intervention NES (GSEA from section 4):
# - **reversed:** significant in both (padj < 0.05), opposite signs;
# - **mimics aging:** significant in both, same sign;
# - **not significant in both:** all others.
#
# **Figure:** each point is a gene set; x = NES with aging, y = NES with the intervention. Green (upper-left, lower-right) = reversed; orange = same direction as aging. The strongest sets of each class are labelled.
# ====================================================================================================

# 5b step 3: per gene set, aging NES vs intervention NES (GSEA of section 4).
# Class: reversed (both padj < 0.05, opposite sign), mimics aging (same sign), not significant in both
options(repr.plot.width = 14, repr.plot.height = 7)

# Function: per pathway, NES of aging next to NES of one intervention, with a class
nes_compare <- function(gs_list, interv) {
  cols <- c("pathway", "NES", "padj")
  m <- merge(as.data.frame(gs_list$old)[, cols], as.data.frame(gs_list[[interv]])[, cols],
             by = "pathway", suffixes = c("_age", "_int"))
  both <- m$padj_age < 0.05 & m$padj_int < 0.05
  both[is.na(both)] <- FALSE
  m$class <- "not significant in both"
  m$class[both & sign(m$NES_age) != sign(m$NES_int)] <- "reversed"
  m$class[both & sign(m$NES_age) == sign(m$NES_int)] <- "mimics aging"
  m$intervention <- c(diet = "aDR vs AL", plasma = "YMP vs PBS")[[interv]]
  m
}

# Function: scatter plot of aging NES vs intervention NES
plot_nes <- function(df, title, prefix, n_labels = 10) {
  df$label <- tolower(gsub("_", " ", sub(prefix, "", df$pathway)))
  lab <- df %>% filter(class != "not significant in both") %>%
    group_by(intervention, class) %>% slice_max(abs(NES_age * NES_int), n = n_labels)
  ggplot(df, aes(x = NES_age, y = NES_int, color = class)) +
    geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
    geom_point(alpha = 0.7) +
    geom_text_repel(data = lab, aes(label = label), size = 2.8, max.overlaps = 30, show.legend = FALSE) +
    scale_color_manual(values = c(reversed = "#1b9e77", `mimics aging` = "#d95f02",
                                  `not significant in both` = "grey80")) +
    facet_wrap(~ intervention) +
    labs(title = title, x = "NES aging (old vs young)", y = "NES intervention", color = NULL) +
    theme_minimal()
}

nes_hallmark <- rbind(nes_compare(gsea_hallmark, "diet"), nes_compare(gsea_hallmark, "plasma"))
nes_gobp     <- rbind(nes_compare(gsea_gobp,     "diet"), nes_compare(gsea_gobp,     "plasma"))

# Number of pathways per class
table(nes_hallmark$intervention, nes_hallmark$class)
table(nes_gobp$intervention,     nes_gobp$class)

plot_nes(nes_hallmark, "Hallmark pathways: aging vs intervention", "HALLMARK_")
plot_nes(nes_gobp,     "GO BP terms: aging vs intervention",       "GOBP_")

# Reversed and aging-mimicking Hallmark pathways
nes_hallmark %>% filter(class != "not significant in both") %>%
  arrange(intervention, class, padj_int) %>%
  select(intervention, class, pathway, NES_age, NES_int, padj_int)

# ====================================================================================================
# #### Step 3b · Reversed pathways at a glance
# - For each reversed set from step 3, the aging NES is shown next to the intervention NES.
# - Hallmark: all reversed pathways, plus those mimicking aging (orange bars); both interventions side by side on a shared pathway axis. No bars in a panel = not significantly changed by that intervention.
# - GO BP: the 10 reversed terms with the lowest intervention padj per intervention.
# - Saved as `reversed_pathways_hallmark.png` and `reversed_pathways_gobp.png`.
#
# **Reading:** grey and blue bars on opposite sides of zero = the intervention pushes the pathway back towards young; longer bars = stronger shift.
# ====================================================================================================

# 5b step 3b: bar figures of reversed (and aging-mimicking) gene sets, aging NES next to intervention NES;
# saved as reversed_pathways_hallmark.png and reversed_pathways_gobp.png

# Function: bar figure of reversed gene sets from nes_compare() output, saved as PNG
# n_top = number of reversed sets shown per intervention (ranked by intervention padj)
# side_by_side = TRUE: interventions as columns with a shared pathway axis
# include_mimics = TRUE: also show gene sets the intervention changes in the same direction as aging
plot_reversed <- function(df, prefix, title, file, n_top = Inf, side_by_side = FALSE, include_mimics = FALSE) {
  classes <- if (include_mimics) c("reversed", "mimics aging") else "reversed"
  r <- df[df$class %in% classes, ]
  if (nrow(r) == 0) { message(title, ": no gene sets to show"); return(invisible(NULL)) }
  r <- do.call(rbind, lapply(split(r, r$intervention), function(d) head(d[order(d$padj_int), ], n_top)))
  r$intervention <- c("aDR vs AL"  = "Dietary restriction (aDR vs AL)",
                      "YMP vs PBS" = "Young plasma (YMP vs PBS)")[r$intervention]

  # readable, wrapped names; pathways ordered by aging NES (same order in both panels)
  r$label <- vapply(tolower(gsub("_", " ", sub(prefix, "", r$pathway))),
                    function(s) paste(strwrap(s, 40), collapse = "\n"), "")
  r$label <- factor(r$label, levels = unique(r$label[order(r$NES_age)]))

  int_lab <- c(reversed = "Intervention: opposite to aging", `mimics aging` = "Intervention: same direction as aging")
  long <- rbind(data.frame(label = r$label, intervention = r$intervention, comparison = "Aging (old vs young)", NES = r$NES_age),
                data.frame(label = r$label, intervention = r$intervention, comparison = int_lab[r$class], NES = r$NES_int))
  long$comparison <- factor(long$comparison, levels = c(rev(int_lab), "Aging (old vs young)"))

  p <- ggplot(long, aes(x = NES, y = label, fill = comparison)) +
    geom_vline(xintercept = 0, color = "grey40") +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = sprintf("%+.2f", NES), hjust = ifelse(NES > 0, -0.2, 1.2)),
              position = position_dodge(width = 0.8), size = 3, color = "grey20") +
    scale_fill_manual(values = setNames(c("grey70", "#3182bd", "#e6550d"), c("Aging (old vs young)", int_lab)),
                      breaks = c("Aging (old vs young)", unname(int_lab)), name = NULL) +
    scale_x_continuous(breaks = seq(-3, 3, 1), expand = expansion(mult = if (side_by_side) 0 else 0.25)) +
    labs(title = title,
         subtitle = if (include_mimics) "Gene sets significant in aging and in the intervention (GSEA padj < 0.05)"
                    else "Significant in aging and in the intervention (GSEA padj < 0.05) with opposite direction",
         x = "Normalized enrichment score (NES)", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top", legend.justification = "left",
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())

  if (side_by_side) {
    # rows where the intervention has no significant effect
    missing <- expand.grid(label = levels(r$label), intervention = unique(r$intervention), stringsAsFactors = FALSE)
    missing <- missing[!paste(missing$label, missing$intervention) %in% paste(r$label, r$intervention), ]
    missing$label <- factor(missing$label, levels = levels(r$label))
    lim <- max(abs(long$NES)) * 1.3   # symmetric axis with room for the value labels
    p <- p + geom_label(data = missing, aes(x = 0, y = label),
                        label = if (include_mimics) "no significant change" else "not reversed", inherit.aes = FALSE,
                        size = 3.2, color = "grey55", fontface = "italic", fill = "white", linewidth = 0) +
      coord_cartesian(xlim = c(-lim, lim)) +
      facet_grid(. ~ intervention) +
      theme(strip.text = element_text(face = "bold", size = 12), panel.spacing.x = unit(2, "lines"))
    fig_w <- 13; fig_h <- 2.2 + 0.55 * nlevels(r$label)
  } else {
    p <- p + facet_grid(intervention ~ ., scales = "free_y", space = "free_y") +
      theme(strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"))
    fig_w <- 10; fig_h <- 2 + 0.55 * nrow(r)
  }
  options(repr.plot.width = fig_w, repr.plot.height = fig_h)
  print(p)
  ggsave(file, p, width = fig_w, height = fig_h, dpi = 300, bg = "white")
  invisible(p)
}

plot_reversed(nes_hallmark, "HALLMARK_", "Hallmark pathways changed by aging and by the interventions",
              "reversed_pathways_hallmark.png", side_by_side = TRUE, include_mimics = TRUE)
plot_reversed(nes_gobp, "GOBP_", "GO BP terms reversed by the interventions (top 10 per intervention)",
              "reversed_pathways_gobp.png", n_top = 10)

# ====================================================================================================
# #### Step 4 · Gene level: which aging genes are reversed?
# Names the genes driving the reversal. Single-gene tests have little power with 3–5 mice, so this level illustrates rather than proves.
# - **Reversed gene:** an aging gene that the intervention changes significantly (padj < 0.05) in the opposite direction.
# - **GO enrichment** of the reversed genes (over-representation, fgsea `fora`). **Background = aging genes of the same direction**, since reversed genes can only come from these.
#
# With small lists, few or no terms are expected at padj < 0.05.
# ====================================================================================================

# 5b step 4: aging genes reversed by an intervention (opposite sign AND intervention padj < 0.05),
# then GO over-representation with background = aging genes of the same direction
age_genes2 <- v2[v2$gene_id %in% unlist(aging_sets), ]
age_genes2$dir <- ifelse(age_genes2$old > 0, "UP with age", "DOWN with age")

reversed2 <- list()
for (i in c("diet", "plasma")) for (d in c("UP with age", "DOWN with age")) {
  g <- age_genes2[age_genes2$dir == d, ]
  p <- g[[paste0("padj_", i)]]
  keep <- sign(g[[i]]) != sign(g$old) & !is.na(p) & p < 0.05
  reversed2[[paste0(i, ", genes ", d)]] <- g[keep, ][order(p[keep]), ]
}

# Number of aging genes and reversed genes per list
data.frame(aging_genes = sapply(names(reversed2), function(n) sum(age_genes2$dir == sub(".*genes ", "", n))),
           reversed    = sapply(reversed2, nrow))

# Top reversed genes per list (log2FC: shrunken)
lapply(names(reversed2), function(n) {
  i <- sub(",.*", "", n)
  head(reversed2[[n]][, c("symbol", "old", i, paste0("padj_", i))], 10)
}) |> setNames(names(reversed2))

# GO enrichment of the reversed genes
# Background = aging genes of the same direction, because reversed genes can only come from those
ora2 <- lapply(names(reversed2), function(n) {
  univ <- unique(na.omit(age_genes2$symbol[age_genes2$dir == sub(".*genes ", "", n)]))
  g    <- unique(na.omit(reversed2[[n]]$symbol))
  if (length(g) < 3) return("fewer than 3 reversed genes, no enrichment test")
  r <- fora(pathways = gobp_sets, genes = g, universe = univ, minSize = 5)
  head(r[order(r$padj), c("pathway", "overlap", "size", "padj")], 10)
})
names(ora2) <- names(reversed2)
ora2

# ====================================================================================================
# #### Step 4b · Gene sets of the reversed genes
# Descriptive, no test, so still informative for short lists.
# - For each reversed gene: its Hallmark and GO BP sets (15–500 genes, the GSEA range). Shown: the 5 most specific (smallest) terms, plus all terms classified as reversed for that intervention in step 3.
# - Outputs: fold-change figure (`aging_genes_reversed_by_interventions.png`; a reversed gene has bars on opposite sides of zero), gene set table (`reversed_genes_gene_sets.png`; full table `reversed_genes_gene_sets.csv`).
#
# **Reading:** gene sets recurring across reversed genes point to restored processes; terms also reversed in step 3 link gene and pathway level.
# ====================================================================================================

# 5b step 4b (descriptive, no test): Hallmark and GO BP sets of each reversed gene;
# outputs a fold-change bar figure and a gene set table (PNG + CSV)

# Function: for every gene symbol, the names of the gene sets that contain it
genes_to_sets <- function(sets) {
  d <- data.frame(symbol = unlist(sets, use.names = FALSE), set = rep(names(sets), lengths(sets)))
  split(d$set, d$symbol)
}
hallmark_of <- genes_to_sets(hallmark_sets)
gobp_of     <- genes_to_sets(gobp_sets)
go_size     <- lengths(gobp_sets)

# Function: readable gene set names (remove prefix, lower case, spaces)
nice <- function(x, prefix) tolower(gsub("_", " ", sub(prefix, "", x)))

# Function: gene set table for one list of reversed genes
# GO BP: only terms with 15-500 genes (same range as GSEA); the 5 most specific (smallest) terms are shown,
# plus all terms that were classified as "reversed" for this intervention in step 3
set_table <- function(n) {
  r <- reversed2[[n]]
  if (nrow(r) == 0) return(NULL)
  i      <- sub(",.*", "", n)
  go_rev <- nes_gobp$pathway[nes_gobp$intervention == c(diet = "aDR vs AL", plasma = "YMP vs PBS")[[i]] &
                             nes_gobp$class == "reversed"]
  do.call(rbind, lapply(seq_len(nrow(r)), function(k) {
    sym <- r$symbol[k]
    hm  <- if (is.na(sym)) character(0) else hallmark_of[[sym]]
    go  <- if (is.na(sym)) character(0) else gobp_of[[sym]]
    go  <- go[go_size[go] >= 15 & go_size[go] <= 500]
    go  <- go[order(go_size[go])]
    data.frame(list                 = n,
               symbol               = sym,
               lfc_aging            = round(r$old[k], 2),
               lfc_intervention     = round(r[[i]][k], 2),
               padj_aging           = signif(r$padj_old[k], 2),
               padj_intervention    = signif(r[[paste0("padj_", i)]][k], 2),
               hallmark             = paste(nice(hm, "HALLMARK_"), collapse = "; "),
               n_go_terms           = length(go),
               go_most_specific     = paste(nice(head(go, 5), "GOBP_"), collapse = "; "),
               go_reversed_in_step3 = paste(nice(intersect(go, go_rev), "GOBP_"), collapse = "; "))
  }))
}

reversed_sets <- do.call(rbind, lapply(names(reversed2), set_table))

# Save as table for the report (opens in Excel)
write.csv(reversed_sets, "reversed_genes_gene_sets.csv", row.names = FALSE)

# ---- Presentable output: fold-change figure + gene set table (both also saved as PNG) ----
library(grid)
library(gridExtra)

has_rev <- !is.null(reversed_sets) && nrow(reversed_sets) > 0
if (!has_rev) {
  message("No aging genes are significantly reversed (padj < 0.05), nothing to show")
} else {
  rs <- reversed_sets
  rs$gene         <- ifelse(is.na(rs$symbol), "(no symbol)", rs$symbol)
  is_diet         <- startsWith(rs$list, "diet")
  rs$intervention <- ifelse(is_diet, "Dietary restriction (vs ad libitum)", "Young plasma (vs PBS)")
  rs$age_dir      <- ifelse(grepl("UP", rs$list), "up with age", "down with age")
  rs$panel        <- ifelse(is_diet, "Dietary restriction vs ad libitum", "Young plasma vs PBS")
  rs$gene_label   <- rs$gene

  # 1) Fold changes: aging bar next to intervention bar; reversal = bars on opposite sides of 0
  lfc_long <- rbind(data.frame(gene = rs$gene_label, panel = rs$panel, comparison = "Aging (old vs young)", lfc = rs$lfc_aging, padj = rs$padj_aging),
                    data.frame(gene = rs$gene_label, panel = rs$panel, comparison = "Intervention (treated vs control)", lfc = rs$lfc_intervention, padj = rs$padj_intervention))
  lfc_long$comparison <- factor(lfc_long$comparison, levels = c("Intervention (treated vs control)", "Aging (old vs young)"))

  fig_h <- 1.6 + 0.8 * nrow(rs) + 0.7 * length(unique(rs$panel))
  p_lfc <- ggplot(lfc_long, aes(x = lfc, y = gene, fill = comparison)) +
    geom_vline(xintercept = 0, color = "grey40") +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = sprintf("%+.2f  (padj = %s)", lfc, formatC(padj, format = "g", digits = 2)),
                  hjust = ifelse(lfc > 0, -0.08, 1.08)),
              position = position_dodge(width = 0.8), size = 3.3, color = "grey20") +
    scale_fill_manual(values = c("Aging (old vs young)" = "grey70", "Intervention (treated vs control)" = "#3182bd"),
                      breaks = c("Aging (old vs young)", "Intervention (treated vs control)"), name = NULL) +
    scale_x_continuous(expand = expansion(mult = 0.8)) +   # room for the labels
    facet_wrap(~ panel, ncol = 1, scales = "free_y", space = "free_y") +   # intervention name above its bars
    labs(title = "Aging genes significantly reversed by an intervention",
         x = "log2 fold change", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top", legend.justification = "left",
          plot.title.position = "plot", legend.location = "plot",
          axis.text.y = element_text(colour = "grey10", size = 15, face = "bold"),
          strip.text = element_text(hjust = 0, face = "bold", size = 12),
          strip.background = element_rect(fill = "grey93", colour = NA),
          panel.spacing = unit(3, "lines"),   # space between the interventions
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
  ggsave("aging_genes_reversed_by_interventions.png", p_lfc, width = 10, height = fig_h, dpi = 150, bg = "white")

  # 2) Gene set table
  wrap <- function(x, width) vapply(x, function(s) paste(strwrap(s, width), collapse = "\n"), "")
  shorten <- function(s, n = 5) {               # first n terms (most specific first) + "(+k more)"
    t <- if (is.na(s) || s == "") character(0) else strsplit(s, "; ", fixed = TRUE)[[1]]
    if (length(t) == 0) return("-")
    paste0(paste(head(t, n), collapse = "; "), if (length(t) > n) sprintf(" (+%d more)", length(t) - n) else "")
  }
  tab <- data.frame(
    "Gene"                     = rs$gene,
    "Intervention"             = rs$intervention,
    "Change with age"          = rs$age_dir,
    "log2FC\naging"            = sprintf("%+.2f", rs$lfc_aging),
    "log2FC\nintervention"     = sprintf("%+.2f", rs$lfc_intervention),
    "padj\nintervention"       = formatC(rs$padj_intervention, format = "g", digits = 2),
    "Hallmark gene sets"       = wrap(ifelse(rs$hallmark == "", "-", rs$hallmark), 30),
    "GO BP terms also reversed\nin step 3 (most specific first)" = wrap(vapply(rs$go_reversed_in_step3, shorten, ""), 45),
    check.names = FALSE)

  tg <- tableGrob(tab, rows = NULL, theme = ttheme_minimal(
    base_size = 10,
    core    = list(fg_params = list(hjust = 0, x = 0.04, vjust = 0.5),
                   bg_params = list(fill = rep(c("white", "grey95"), length.out = nrow(tab)))),
    colhead = list(fg_params = list(hjust = 0, x = 0.04, fontface = "bold"),
                   bg_params = list(fill = "grey85"))))
  tab_w <- convertWidth(sum(tg$widths),  "in", valueOnly = TRUE) + 0.3
  tab_h <- convertHeight(sum(tg$heights), "in", valueOnly = TRUE) + 0.3
  png("reversed_genes_gene_sets.png", width = tab_w, height = tab_h, units = "in", res = 150, bg = "white")
  grid.draw(tg); invisible(dev.off())
}

# Draw the figure and the table as separate steps: the notebook draws a plot only after the
# whole if-block has finished, so each needs its own block to be shown at its own size
if (has_rev) {
  options(repr.plot.width = 10, repr.plot.height = fig_h)
  print(p_lfc)
}
if (has_rev) {
  options(repr.plot.width = tab_w, repr.plot.height = tab_h)
  grid.newpage(); grid.draw(tg)
}

# ====================================================================================================
# #### Control 1 · Specificity for aging genes
# Fold changes from two datasets can agree for reasons unrelated to aging (technical biases, general stress responses). A rejuvenating effect should be stronger for aging genes than for random genes.
# - Spearman correlation of aging vs intervention shrunken log2FC over the aging genes, compared with 1,000 random sets of non-aging genes of equal size, matched by expression level (deciles of mean normalized counts).
# - Empirical two-sided p-value: share of random sets deviating at least as much from the random mean; Benjamini–Hochberg across the two interventions.
#
# **Reading:** observed correlation (red line) outside the random distribution (grey) = specific to aging genes; negative = rejuvenation, positive = mimics aging.
# ====================================================================================================

# 5b control 1: is the agreement with aging specific to aging genes?
# Score = Spearman rho between aging and intervention shrunken fold changes over the aging genes.
# Null = the same score for 1000 random gene sets of the same size, matched for expression level.
options(repr.plot.width = 10, repr.plot.height = 4)

v2$expr_bin <- cut(rank(v2$baseMean), 10, labels = FALSE)   # expression deciles
sig_ids   <- unlist(aging_sets)
pool      <- v2[!v2$gene_id %in% sig_ids, ]
pool_bins <- split(seq_len(nrow(pool)), pool$expr_bin)
need      <- table(v2$expr_bin[v2$gene_id %in% sig_ids])     # genes needed per decile

# Function: observed rho over the aging genes vs n_perm random expression-matched sets; empirical two-sided p
null_test <- function(interv, n_perm = 1000) {
  sig <- v2[v2$gene_id %in% sig_ids, ]
  obs <- cor(sig$old, sig[[interv]], method = "spearman")
  set.seed(42)
  null <- replicate(n_perm, {
    idx <- unlist(lapply(names(need), function(b) sample(pool_bins[[b]], need[[b]])))
    cor(pool$old[idx], pool[[interv]][idx], method = "spearman")
  })
  list(summary = data.frame(rho_aging_genes = round(obs, 3),
                            rho_random_mean = round(mean(null), 3),
                            rho_random_2.5  = round(quantile(null, 0.025), 3),
                            rho_random_97.5 = round(quantile(null, 0.975), 3),
                            p_empirical     = (sum(abs(null - mean(null)) >= abs(obs - mean(null))) + 1) / (n_perm + 1)),
       null = null, obs = obs)
}

null_res <- list(diet = null_test("diet"), plasma = null_test("plasma"))
null_tab <- do.call(rbind, lapply(null_res, `[[`, "summary"))
null_tab$padj <- p.adjust(null_tab$p_empirical, method = "BH")   # adjusted across the two interventions
null_tab$p_empirical <- NULL
null_tab

null_df <- do.call(rbind, lapply(names(null_res), function(n) data.frame(intervention = n, rho = null_res[[n]]$null)))
obs_df  <- data.frame(intervention = names(null_res), rho = sapply(null_res, `[[`, "obs"))
ggplot(null_df, aes(x = rho)) +
  geom_histogram(bins = 40, fill = "grey70") +
  geom_vline(data = obs_df, aes(xintercept = rho), color = "#b2182b", linewidth = 1) +
  facet_wrap(~ intervention) +
  labs(title = "Aging genes vs 1000 random gene sets",
       subtitle = "Red line: observed correlation over aging genes; grey: random, expression-matched sets",
       x = "Spearman rho (aging vs intervention log2FC)", y = "Count") +
  theme_minimal()

# ====================================================================================================
# #### Control 2 · Robustness to the aging-signature definition
# The camera test of step 2, repeated with alternative aging signatures:
# - **strict:** padj < 0.05 and |shrunken log2FC| > 0.5;
# - **18 vs 3** and **21 vs 3 months** separately;
# - **same sex:** female aging genes for aDR, male for YMP (2 young vs 4 old per sex, so smaller, less precise signatures).
#
# **Reading:** robust if the direction (ideally also significance) is consistent across signatures.
# ====================================================================================================

# 5b control 2: the camera test of step 2, repeated with alternative aging signatures
# (strict, 18 vs 3 mo, 21 vs 3 mo, same sex)

# Function: aging-up / aging-down sets from a DE result (padj < 0.05, common genes only)
make_sets <- function(res, min_lfc = 0) {
  r <- res[res$gene_id %in% common_genes & !is.na(res$padj) & res$padj < 0.05 &
           abs(res$log2FoldChange) > min_lfc, ]
  list(aging_up   = r$gene_id[r$log2FoldChange > 0],
       aging_down = r$gene_id[r$log2FoldChange < 0])
}

# Sex-specific aging signatures: diet mice are female, plasma mice are male
fit_aging_sex <- function(sx) {
  keep <- cd_age$sex == sx
  dds <- DESeqDataSetFromMatrix(cts_age[, keep], colData = cd_age[keep, ], design = ~ age_bin)
  get_res(DESeq(dds, quiet = TRUE), c("age_bin", "old", "young"))
}
de_old_female <- cached("de_old_female", fit_aging_sex("female"))
de_old_male   <- cached("de_old_male",   fit_aging_sex("male"))

strict <- v2[which(v2$padj_old < 0.05 & abs(v2$old) > 0.5), ]
signatures <- list(
  "main (old vs young)"           = aging_sets,
  "strict (|shrunken LFC| > 0.5)" = list(aging_up   = strict$gene_id[strict$old > 0],
                                         aging_down = strict$gene_id[strict$old < 0]),
  "18 vs 3 mo only"               = make_sets(de_age18),
  "21 vs 3 mo only"               = make_sets(de_age21)
)

robust <- do.call(rbind, c(
  lapply(names(signatures), function(s) rbind(
    cbind(signature = s, intervention = "aDR vs AL",  run_camera(dds_diet,   signatures[[s]])),
    cbind(signature = s, intervention = "YMP vs PBS", run_camera(dds_plasma, signatures[[s]])))),
  list(cbind(signature = "same sex (female)", intervention = "aDR vs AL",  run_camera(dds_diet,   make_sets(de_old_female))),
       cbind(signature = "same sex (male)",   intervention = "YMP vs PBS", run_camera(dds_plasma, make_sets(de_old_male))))
))
rownames(robust) <- NULL
robust[order(robust$intervention, robust$set), ]

# ====================================================================================================
# #### Control 3 · Cell-type composition
# A difference in tissue composition (e.g. more white matter from dissection) can mimic changes in myelin or oligodendrocyte genes without any change in regulation.
# - Per sample and cell type: mean z-score (within dataset) of log TPM of marker genes: microglia (*P2ry12, Tmem119, Cx3cr1, Hexb, Csf1r*), oligodendrocytes (*Mbp, Plp1, Mog, Mobp, Mag*), astrocytes (*Aqp4, Aldh1l1, Gja1, Slc1a3*), granule neurons (*Gabra6, Neurod1, Cbln3, Rbfox3*), Purkinje cells (*Pcp2, Calb1, Car8*).
# - Treated vs control: Welch's t-test; Benjamini–Hochberg across cell types and interventions.
#
# **Reading:** a clear difference means effects on that cell type's genes may partly reflect composition. Microglial scores can also rise with activation, not only cell number.
# ====================================================================================================

# 5b control 3: do cell-type proportions differ between treated and control mice?
# Score per sample = mean z-score (within its dataset) of log TPM of a few marker genes per cell type
options(repr.plot.width = 12, repr.plot.height = 9)

markers <- list(
  "Microglia"        = c("P2ry12", "Tmem119", "Cx3cr1", "Hexb", "Csf1r"),
  "Oligodendrocytes" = c("Mbp", "Plp1", "Mog", "Mobp", "Mag"),
  "Astrocytes"       = c("Aqp4", "Aldh1l1", "Gja1", "Slc1a3"),
  "Granule neurons"  = c("Gabra6", "Neurod1", "Cbln3", "Rbfox3"),
  "Purkinje cells"   = c("Pcp2", "Calb1", "Car8")
)

# Function: one gene_id per marker symbol (the most highly expressed one)
marker_ids <- function(symbols) {
  sapply(symbols, function(s) {
    ids <- meta_genes$gene_id[which(meta_genes$symbol == s)]
    ids[which.max(rowMeans(expr[ids, , drop = FALSE]))]
  })
}

# Score per sample and cell type: mean z-score of its markers' log(TPM + 1), within each dataset
ct_scores <- do.call(rbind, lapply(unique(meta$dataset), function(ds) {
  s <- meta$SRR[meta$dataset == ds]
  do.call(rbind, lapply(names(markers), function(ct) {
    x <- log1p(expr[marker_ids(markers[[ct]]), s, drop = FALSE])
    z <- t(scale(t(x)))                                  # z-score per gene within this dataset
    data.frame(SRR = s, cell_type = ct, score = colMeans(z, na.rm = TRUE))
  }))
}))
ct_scores <- merge(ct_scores, meta[, c("SRR", "group", "cohort")], by = "SRR")

ggplot(ct_scores, aes(x = group, y = score, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1.5) +
  scale_fill_manual(values = group_cols, guide = "none") +
  facet_grid(cell_type ~ cohort, scales = "free_x") +
  labs(title = "Cell-type marker scores per sample", x = NULL, y = "Mean marker z-score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# t-test per cell type: treated vs control within each intervention cohort
ct_tests <- do.call(rbind, lapply(list(c("AL", "aDR"), c("PBS", "YMP")), function(g) {
  do.call(rbind, lapply(names(markers), function(ct) {
    d <- ct_scores[ct_scores$cell_type == ct, ]
    tt <- t.test(d$score[d$group == g[2]], d$score[d$group == g[1]])
    data.frame(comparison = paste(g[2], "vs", g[1]), cell_type = ct,
               difference = round(unname(tt$estimate[1] - tt$estimate[2]), 2), p = tt$p.value)
  }))
}))
ct_tests$padj <- signif(p.adjust(ct_tests$p, method = "BH"), 2)   # adjusted across all 10 tests
ct_tests$p <- NULL
ct_tests

# ====================================================================================================
# ## 6 · Conclusions
# - **Aging:** more inflammation and immune activation (microglia, complement, interferon response, MHC class I); less myelin/oligodendrocyte and energy-metabolism gene expression. Mostly established by 18 months.
# - **Young plasma:** few individually significant genes, but a broad reversal that is specific to aging genes and strongest for inflammatory microglial genes. Power is limited (n = 3).
# - **Dietary restriction:** a robust but distinct effect, dominated by circadian, glucocorticoid-responsive and metabolic changes. It attenuates the strongest inflammatory genes and restores oxidative phosphorylation, but reinforces the decline of myelin genes, so it is not a general rejuvenating intervention in the cerebellum.
#
# **Limitations**
# - Aging and intervention samples come from separate datasets, so only fold changes can be compared.
# - Small groups (3–5 mice).
# - aDR effects may be confounded by feeding time (circadian genes).
# - Dissection differences may affect tissue composition in the plasma cohort.
# ====================================================================================================

