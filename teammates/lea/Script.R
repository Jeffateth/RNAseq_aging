# Lea's analysis, refactored to the team's common setup. Exported from teammates/lea/Script.ipynb (R); text cells appear as comments.
# Run from teammates/lea/ with the repo's pixi environment.

# ====================================================================================================
# # Aging and rejuvenating interventions in the mouse cerebellum
#
# **Objective:** describe transcriptional aging in the mouse cerebellum, describe the effects of acute dietary restriction (aDR) and young mouse plasma (YMP), and explore whether each intervention opposes the aging changes.
#
# - **Reference:** Hahn et al. (2023). **Data:** RSEM quantification (STAR alignments, GENCODE M39) of 27 runs, BioProjects PRJNA875066 (aging) and PRJNA946653 (interventions), read from the repo's `for_teammates/` folder.
# - **Setup:** the team's common DESeq2 setup (agreed 20 November 2026), summarised in Section 3.
# - **Run:** from this folder with the repo's pixi environment (R kernel). Outputs are written to `output/`.
#
# | Section | Question |
# |---|---|
# | 1 · Import and quality control | Which samples go into the analysis? |
# | 2 · Exploratory analysis | Do samples group by experiment and condition? |
# | 3 · Differential expression | Which genes change with age, aDR and YMP? |
# | 4 · Pathways (Hallmark GSEA) | Which biological processes shift? |
# | 5 · Aging vs interventions | Do aDR or YMP move aging genes in the opposite direction? |
# ====================================================================================================

# Packages, figure style and output folder
suppressPackageStartupMessages({
  library(DESeq2); library(apeglm); library(fgsea); library(msigdbr); library(ggrepel); library(cowplot)
  library(tidyverse)   # loaded last, so dplyr's count/select/filter are not masked
})
options(jupyter.rich_display = FALSE, width = 140)
theme_set(theme_minimal(base_size = 26) + theme(panel.grid.minor = element_blank(), plot.title = element_text(size = rel(1)),
                                              strip.text = element_text(size = rel(1)), legend.position = "bottom"))
DATA <- "../../for_teammates"   # shared RSEM tables (in git)
OUT  <- "output"; dir.create(OUT, showWarnings = FALSE)
ALPHA <- 0.05                   # DEG and gene set cutoff (padj), also DESeq2's results(alpha = )
COL <- c(aging = "#2a78d6", aDR = "#eb6834", YMP = "#1baf7a", up = "#e34948", down = "#2a78d6", ns = "#cfcfcf")
save_png <- function(p, name, w, h) ggsave(file.path(OUT, name), p, width = w, height = h, dpi = 150, bg = "white")

# ====================================================================================================
# ## 1 · Import and quality control
# Three RSEM tables per gene and sample: **expected counts** (reads), **TPM**, and **effective length** (the gene length RSEM used, per sample). DESeq2 receives them the way `tximport` would ("tximport workflow"): counts plus effective lengths, so that its normalisation corrects for gene length differences between samples as well as for sequencing depth.
# ====================================================================================================

# Load the three RSEM tables (genes x samples) and check that they match
rd <- function(f) as.matrix(read.csv(file.path(DATA, f), row.names = 1, check.names = FALSE))
counts <- rd("rsem_expected_counts.csv.gz"); tpm <- rd("rsem_tpm.csv.gz"); len <- rd("rsem_effective_length.csv.gz")
stopifnot(identical(dimnames(counts), dimnames(tpm)), identical(dimnames(counts), dimnames(len)))
genes <- read.csv(file.path(DATA, "genes_M39.csv"))
sym <- setNames(genes$symbol, genes$gene_id)

# Effective lengths must be present for every gene and sample (no other normalisation is substituted)
stopifnot(!anyNA(len))
zero <- len <= 0
cat(nrow(counts), "genes x", ncol(counts), "samples\n")
cat("Effective length 0:", sum(zero), "entries; their counts:", sum(counts[zero]), "(genes without reads in that sample)\n")
len[zero] <- 1   # DESeq2 needs lengths > 0; these entries have no reads, so the value has no effect
txi <- list(abundance = tpm, counts = counts, length = len, countsFromAbundance = "no")

# Sample metadata: one row per run, in the table's column order
meta <- read.csv(file.path(DATA, "samples.csv")) %>%
  mutate(experiment = case_when(cohort == "Aging" ~ "aging", group %in% c("AL", "aDR") ~ "aDR", TRUE ~ "YMP"),
         group = factor(group, levels = c("aging 3mo", "aging 18mo", "aging 21mo", "AL", "aDR", "PBS", "YMP")),
         sex = factor(sex))
rownames(meta) <- meta$SRR
meta <- meta[colnames(counts), ]
stopifnot(identical(rownames(meta), colnames(counts)))

# ====================================================================================================
# ### Sample exclusion
# **SRR21355389** (aging, 21 months, female) is excluded: its haemoglobin genes are far higher than in every other aging sample (table below), indicating blood contamination of the tissue. **SRR23914572** (AL) is kept, as agreed: its unusual genes are choroid plexus genes (a neighbouring tissue), not blood.
# ====================================================================================================

# Haemoglobin TPM in the aging samples: SRR21355389 vs the other 11
hb <- names(sym)[sym %in% c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt")]
aging_srr <- meta$SRR[meta$experiment == "aging"]
others <- setdiff(aging_srr, "SRR21355389")
tibble(gene = sym[hb],
       SRR21355389 = round(tpm[hb, "SRR21355389"]),
       others_median = round(apply(tpm[hb, others], 1, median)),
       others_max = round(apply(tpm[hb, others], 1, max)),
       x_vs_max = round(tpm[hb, "SRR21355389"] / apply(tpm[hb, others], 1, max), 1))

# Exclude SRR21355389 everywhere; sample counts per group after exclusion
EXCLUDE <- "SRR21355389"
meta <- meta[!rownames(meta) %in% EXCLUDE, ]
txi[c("abundance", "counts", "length")] <- lapply(txi[c("abundance", "counts", "length")], \(m) m[, rownames(meta)])
meta %>% count(experiment, group, sex) %>% pivot_wider(names_from = sex, values_from = n, values_fill = 0) %>%
  mutate(total = female + male)

# ====================================================================================================
# ## 2 · Exploratory analysis
# - **Genes:** at least 10 reads in at least 3 samples (3 = the smallest group).
# - **Transformation:** variance-stabilising transformation (`vst`, blind to the groups): log-like counts, corrected for sequencing depth and gene length, in which low and high counts vary similarly.
# - **PCA:** on the 500 most variable genes.
# ====================================================================================================

# VST of all 26 samples, then PCA (top 500 genes); colour = experiment, shape = condition
vst_of <- function(srr, min_n = 3) {
  x <- txi; x[c("abundance", "counts", "length")] <- lapply(x[c("abundance", "counts", "length")], \(m) m[, srr])
  d <- suppressMessages(DESeqDataSetFromTximport(x, meta[srr, ], ~ 1))
  vst(d[rowSums(counts(d) >= 10) >= min_n, ], blind = TRUE)
}
pca_df <- function(v) {
  p <- plotPCA(v, intgroup = "group", returnData = TRUE, ntop = 500)
  list(d = cbind(p, meta[rownames(p), c("experiment", "sex")]), pv = round(100 * attr(p, "percentVar")))
}
v_all <- vst_of(rownames(meta))
p <- pca_df(v_all)
p_all <- ggplot(p$d, aes(PC1, PC2, color = experiment, shape = experiment)) + geom_point(size = 7, alpha = 0.85) +
  scale_color_manual(values = COL) + labs(x = paste0("PC1 (", p$pv[1], "%)"), y = paste0("PC2 (", p$pv[2], "%)"),
                                          color = NULL, shape = NULL, title = "All samples")
options(repr.plot.width = 10, repr.plot.height = 8.5); p_all
save_png(p_all, "pca_all.png", 10, 8.5)

# PCA within each experiment (each with its own VST and top 500 genes); colour = condition
within <- lapply(c("aging", "aDR", "YMP"), \(e) {
  p <- pca_df(vst_of(rownames(meta)[meta$experiment == e]))
  p$d$condition <- droplevels(p$d$group)
  ggplot(p$d, aes(PC1, PC2, color = condition, shape = sex)) + geom_point(size = 6) +
    scale_color_brewer(palette = "Dark2") + labs(title = e, color = NULL, shape = NULL,
      x = paste0("PC1 (", p$pv[1], "%)"), y = paste0("PC2 (", p$pv[2], "%)")) +
    guides(color = guide_legend(nrow = 3), shape = guide_legend(nrow = 2)) + theme(legend.text = element_text(size = 18))
})
p_within <- plot_grid(plotlist = within, nrow = 1)
options(repr.plot.width = 24, repr.plot.height = 9); p_within
save_png(p_within, "pca_within.png", 24, 9)

# **Reading the PCA**
# - **All samples:** PC1 (41 %) separates the aging cohort from the two intervention cohorts, PC2 (13 %) separates aDR from YMP. The largest differences lie between experiments (separate studies, different sexes and ages), so each comparison is analysed within its own cohort only.
# - **Within the aging cohort:** PC1 (28 %) mainly separates males from females, more than age does; this is why the aging model includes a sex term.
# - **Within aDR and YMP:** treated and control mice partly separate but overlap; with 3–5 mice per group, the treatment effects are modest.

# ====================================================================================================
# ## 3 · Differential expression (DESeq2)
# | Choice | Setting (team setup) |
# |---|---|
# | Input | RSEM expected counts + effective lengths (`DESeqDataSetFromTximport`) |
# | Aging | old (18 + 21 months, n = 7) vs young (3 months, n = 4), design `~ sex + age_bin` |
# | aDR | aDR (5) vs AL (4), design `~ condition` (all female) |
# | YMP | YMP (3) vs PBS (3), design `~ condition` (all male) |
# | Gene filter | ≥ 10 reads in at least *n* samples, *n* = smallest group (4, 4, 3) |
# | Test | Wald, `results(alpha = 0.05)` |
# | DEG | padj < 0.05, no fold-change cutoff |
# | Fold changes in plots and comparisons | apeglm-shrunken (`lfcShrink`); p-values are not shrunken |
#
# **Tested genes** (team definition) = genes that pass the gene filter and so get a Wald statistic, fold change and p-value. Some of them get padj = NA from DESeq2's independent filtering (too few reads to reach padj < 0.05) or a count outlier; they cannot be DEGs but still enter the GSEA ranking.
# ====================================================================================================

# Fit one DESeq2 model: samples, design, the tested factor (with the control level first), min_n for the filter
fit_de <- function(srr, design, factor, levels, min_n) {
  cd <- meta[srr, ]; cd[[factor]] <- factor(cd[[factor]], levels = levels)
  x <- txi; x[c("abundance", "counts", "length")] <- lapply(x[c("abundance", "counts", "length")], \(m) m[, srr])
  d <- suppressMessages(DESeqDataSetFromTximport(x, cd, design))
  d <- d[rowSums(counts(d) >= 10) >= min_n, ]
  DESeq(d, quiet = TRUE)
}
# Results table for the tested factor: Wald test (alpha = 0.05) + apeglm-shrunken log2FC + gene symbol, sorted by padj
de_table <- function(dds) {
  coef <- resultsNames(dds)[length(resultsNames(dds))]   # last coefficient = the tested factor
  r <- as.data.frame(results(dds, name = coef, alpha = ALPHA))
  s <- lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE)
  tibble(gene_id = rownames(r), symbol = sym[rownames(r)], baseMean = r$baseMean, log2FoldChange = r$log2FoldChange,
         lfc_shrunk = s$log2FoldChange, stat = r$stat, pvalue = r$pvalue, padj = r$padj) %>% arrange(padj)
}

meta$age_bin   <- ifelse(meta$age_months == 3 & meta$experiment == "aging", "young", "old")
meta$condition <- as.character(meta$group)
srr <- split(rownames(meta), meta$experiment)
dds <- list(
  aging = fit_de(srr$aging, ~ sex + age_bin, "age_bin",   c("young", "old"), min_n = 4),
  aDR   = fit_de(srr$aDR,   ~ condition,     "condition", c("AL", "aDR"),    min_n = 4),
  YMP   = fit_de(srr$YMP,   ~ condition,     "condition", c("PBS", "YMP"),   min_n = 3))
res <- lapply(dds, de_table)
sapply(dds, \(d) tail(resultsNames(d), 1))   # tested coefficient per comparison

# DEG summary (reference from the team setup: aging 175, aDR 118, YMP 50)
deg_summary <- imap_dfr(res, \(r, n) tibble(comparison = n, tested = sum(!is.na(r$stat)), with_padj = sum(!is.na(r$padj)),
  DEGs = sum(r$padj < ALPHA, na.rm = TRUE), up = sum(r$padj < ALPHA & r$log2FoldChange > 0, na.rm = TRUE),
  down = sum(r$padj < ALPHA & r$log2FoldChange < 0, na.rm = TRUE))) %>%
  mutate(reference = c(175, 118, 50))
deg_summary

# Volcano plots: x = shrunken log2FC, colour = DEG (padj < 0.05) up / down; the 8 most significant DEGs are named
vd <- bind_rows(res, .id = "comparison") %>% filter(!is.na(padj)) %>%   # genes with padj only (y axis)
  mutate(comparison = factor(comparison, names(res)),
         call = case_when(padj < ALPHA & log2FoldChange > 0 ~ "up", padj < ALPHA ~ "down", TRUE ~ "ns"),
         title = paste0(comparison, ": ", ave(call == "up", comparison, FUN = sum), " up, ", ave(call == "down", comparison, FUN = sum), " down"),
         title = fct_reorder(title, as.integer(comparison)))
top <- vd %>% filter(call != "ns") %>% group_by(comparison) %>% slice_min(padj, n = 8, with_ties = FALSE)
p_volc <- ggplot(vd, aes(lfc_shrunk, -log10(padj), color = call)) + geom_point(size = 1.5) +
  geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", color = "grey50") +
  geom_text_repel(data = top, aes(label = symbol), color = "black", size = 6, seed = 42, max.overlaps = 30) +
  scale_color_manual(values = COL, breaks = c("up", "down"), labels = c("up", "down")) +
  facet_wrap(~ title, scales = "free") + labs(x = "log2 fold change (shrunken)", y = "-log10 padj", color = NULL) +
  guides(color = guide_legend(override.aes = list(size = 6)))
options(repr.plot.width = 24, repr.plot.height = 9); p_volc
save_png(p_volc, "volcano.png", 24, 9)

# Top 10 DEGs per comparison (shrunken log2FC; positive = higher in old / treated mice)
bind_rows(res, .id = "comparison") %>% filter(padj < ALPHA) %>% group_by(comparison) %>%
  slice_min(padj, n = 10, with_ties = FALSE) %>%
  summarise(top10 = paste0(symbol, " (", sprintf("%+.1f", lfc_shrunk), ")", collapse = ", ")) %>%
  mutate(comparison = factor(comparison, names(res))) %>% arrange(comparison)

# Heatmap of the 30 most significant aging DEGs: VST of the aging samples, z-score per gene; columns = mice by age and sex
v_age <- vst_of(srr$aging, min_n = 4)
top_age <- res$aging %>% filter(padj < ALPHA) %>% slice_min(padj, n = 30, with_ties = FALSE)
z <- t(scale(t(assay(v_age)[top_age$gene_id, ])))
ord_s <- rownames(meta[srr$aging, ] %>% arrange(age_months, sex))
ord_g <- top_age$gene_id[order(top_age$log2FoldChange)]
hd <- as.data.frame(as.table(z)) %>% setNames(c("gene_id", "SRR", "z")) %>%
  mutate(gene = factor(sym[as.character(gene_id)], unique(sym[ord_g])),
         mouse = factor(paste0(meta[as.character(SRR), "age_months"], "mo ", substr(meta[as.character(SRR), "sex"], 1, 1)),
                        levels = unique(paste0(meta[ord_s, "age_months"], "mo ", substr(meta[ord_s, "sex"], 1, 1)))),
         SRR = factor(SRR, ord_s), z = pmax(pmin(z, 2.5), -2.5))
p_heat <- ggplot(hd, aes(SRR, gene, fill = z)) + geom_tile(color = "white") +
  scale_fill_gradient2(low = COL[["down"]], mid = "#f4f4f4", high = COL[["up"]], limits = c(-2.5, 2.5)) +
  scale_x_discrete(labels = setNames(as.character(hd$mouse), hd$SRR)) +
  labs(x = NULL, y = NULL, fill = "z-score", title = "Top 30 aging DEGs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 20), axis.text.y = element_text(face = "italic", size = 18),
        panel.grid = element_blank(), legend.position = "right", legend.key.height = unit(2.5, "lines"))
options(repr.plot.width = 11, repr.plot.height = 13); p_heat
save_png(p_heat, "heatmap_aging.png", 11, 13)

# Save the three DE tables (all tested genes; padj NA = removed by independent filtering or a count outlier)
iwalk(res, \(r, n) write_csv(r, file.path(OUT, paste0("de_", n, ".csv"))))
write_csv(deg_summary, file.path(OUT, "deg_summary.csv"))

# **Reading the differential expression**
# - **DEG counts** reproduce the team reference exactly: aging 175, aDR 118, YMP 50.
# - **Aging (138 up, 37 down):** the top genes form a neuroinflammation signature: complement (*C4b, C1qa, C1qb, C1qc*), activated microglia (*Clec7a, Cst7, Trem2, Tyrobp, Cd68, Lgals3, Gpnmb*), MHC class I (*H2-K1, H2-Q6, B2m*), lysosomal cathepsins (*Ctss, Ctsd, Ctsz*) and the astrocyte gene *Gfap*. In the heatmap, all 7 old mice (18 and 21 months alike) are above the 4 young mice for these genes.
# - **aDR (46 up, 72 down):** mostly small fold changes. Top genes include glucocorticoid- and circadian-responsive genes (*Tsc22d3, Hlf, Ciart*) and immediate early genes (*Nr4a2, Nr4a3*), not the aging inflammation genes. *Cdr1* (log2FC +12) has 0 reads in all 4 AL mice and 15–555 in the aDR mice (it is also detected in some PBS / YMP mice): an on/off pattern to interpret with caution.
# - **YMP (22 up, 28 down):** few genes; several large fold changes rest on 3 vs 3 mice (*Ccn3, Slc47a1, Gm12895*). The most significant gene, *Ptgds*, and *Ogn* and *Slc47a1* are typical of the meninges, so part of this signal may reflect how much meninges was included in each dissection rather than a plasma effect.

# ====================================================================================================
# ## 4 · Pathways: Hallmark GSEA
# - **Gene sets:** mouse-native MSigDB Hallmark (`msigdbr`, `db_species = "MM"`, collection MH; 50 sets).
# - **Ranking:** all tested genes, ordered by the DESeq2 Wald statistic (log2FC / its standard error); one value per gene symbol.
# - **Test:** `fgsea` (sets of 15–500 genes, seed 42). **NES** > 0 = set higher in old / treated mice; significant = padj < 0.05.
# ====================================================================================================

# Hallmark GSEA for each comparison
hallmark <- msigdbr(species = "Mus musculus", db_species = "MM", collection = "MH")
hallmark <- split(hallmark$gene_symbol, hallmark$gs_name)
cat("msigdbr", as.character(packageVersion("msigdbr")), "|", length(hallmark), "Hallmark sets\n")
run_gsea <- function(r) {
  r <- r %>% filter(!is.na(stat), !is.na(symbol), symbol != "") %>% group_by(symbol) %>%
    slice_max(abs(stat), n = 1, with_ties = FALSE) %>% ungroup()
  set.seed(42)
  suppressWarnings(fgsea(hallmark, setNames(r$stat, r$symbol), minSize = 15, maxSize = 500)) %>%
    as_tibble() %>% select(pathway, NES, pval, padj, size) %>% arrange(padj)
}
gsea <- bind_rows(lapply(res, run_gsea), .id = "comparison")
write_csv(gsea, file.path(OUT, "gsea_hallmark.csv"))
gsea %>% group_by(comparison = factor(comparison, names(res))) %>%
  summarise(significant = sum(padj < ALPHA), up = sum(padj < ALPHA & NES > 0), down = sum(padj < ALPHA & NES < 0))

# One figure: NES of every Hallmark set significant in at least one comparison (* = padj < 0.05), ordered by aging NES
nice <- \(p) str_to_sentence(gsub("_", " ", sub("^HALLMARK_", "", p))) %>%
  str_replace_all(c("^Tnfa signaling via nfkb" = "TNFa signaling via NF-kB", "^Il6 jak stat3" = "IL6-JAK-STAT3", "^Il2 stat5" = "IL2-STAT5",
                    "^Mtorc1" = "mTORC1", "^Myc" = "MYC", "^E2f" = "E2F", "^G2m" = "G2M", "^Kras" = "KRAS", "^P53" = "p53", "^Uv" = "UV",
                    "^Tgf beta" = "TGF-beta", "^Pi3k akt mtor" = "PI3K-AKT-mTOR", "^Dna" = "DNA", "^Wnt beta catenin" = "WNT-beta-catenin"))
keep <- unique(gsea$pathway[gsea$padj < ALPHA])
ord  <- gsea %>% filter(comparison == "aging") %>% arrange(NES) %>% pull(pathway)
gd <- gsea %>% filter(pathway %in% keep) %>%
  mutate(set = factor(nice(pathway), nice(ord[ord %in% keep])), comparison = factor(comparison, names(res)))
p_gsea <- ggplot(gd, aes(comparison, set, fill = NES)) + geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = ifelse(padj < ALPHA, "*", "")), size = 12, vjust = 0.75) +
  scale_fill_gradient2(low = COL[["down"]], mid = "#f4f4f4", high = COL[["up"]]) + scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, title = "Hallmark gene sets (* padj < 0.05)") +
  theme(panel.grid = element_blank(), legend.position = "right", legend.key.height = unit(2.5, "lines"), axis.text.y = element_text(size = 20))
options(repr.plot.width = 13, repr.plot.height = 15); p_gsea
save_png(p_gsea, "gsea_hallmark.png", 13, 15)

# Hallmark sets significant in aging: how many are also significant in each intervention, in the opposite or the same direction
aging_sets <- gsea %>% filter(comparison == "aging", padj < ALPHA) %>% select(pathway, NES_aging = NES)
gsea %>% filter(comparison != "aging", padj < ALPHA) %>% inner_join(aging_sets, by = "pathway") %>%
  group_by(intervention = comparison) %>%
  summarise(aging_sets = nrow(aging_sets), significant_in_both = n(),
            opposite = sum(sign(NES) != sign(NES_aging)), same = sum(sign(NES) == sign(NES_aging)),
            opposite_sets = paste(nice(pathway[sign(NES) != sign(NES_aging)]), collapse = ", "),
            same_sets = paste(nice(pathway[sign(NES) == sign(NES_aging)]), collapse = ", "))

# **Reading the pathways**
# - **Aging (15 sets):** up: interferon α and γ responses, allograft rejection, complement, inflammatory response, coagulation, KRAS signalling; down: oxidative phosphorylation, mTORC1, fatty acid metabolism, MYC targets, TNFα/NF-κB, glycolysis, hypoxia, androgen response.
# - **aDR (16 sets):** all 6 aging sets that are also significant in aDR move in the opposite direction (complement, coagulation and KRAS down; oxidative phosphorylation, MYC targets and TNFα/NF-κB up).
# - **YMP (24 sets, 22 up):** 9 of the 11 shared sets move opposite to aging (interferon responses down; several metabolic sets up), 2 in the same direction (complement, coagulation). Because almost all YMP sets go up, the 'opposite' calls are mostly the sets that go down with age; this is partly a general upward shift, not a specific reversal.
# - Gene sets share genes, so these counts are not independent pieces of evidence.

# ====================================================================================================
# ## 5 · Aging vs interventions
# The cohorts are separate experiments, so only fold changes are compared, never expression values.
#
# 1. **Aging genes:** aging DEGs (padj < 0.05).
# 2. Of these, the genes with a shrunken fold change in the intervention model.
# 3. Per intervention: aging log2FC (x) vs intervention log2FC (y). Points in the upper-left and lower-right quadrants move **opposite** to aging.
# 4. **Spearman ρ** (negative = opposite to aging overall) and **% opposite sign**, as descriptive numbers. Genes that are also DEGs of the intervention are highlighted.
#
# Not used: a test of whether the % opposite differs from 50 % (aging genes are correlated, not independent, so such a p-value would overstate the evidence).
# ====================================================================================================

# Compare each intervention with aging over the aging DEGs
aging_deg <- res$aging %>% filter(padj < ALPHA) %>% select(gene_id, symbol, aging = lfc_shrunk)
cmp <- bind_rows(lapply(c("aDR", "YMP"), \(n) aging_deg %>%
  inner_join(res[[n]] %>% select(gene_id, lfc_int = lfc_shrunk, padj_int = padj), by = "gene_id") %>%
  filter(!is.na(lfc_int)) %>%
  mutate(panel = n, opposite = sign(aging) != sign(lfc_int), both_DEG = !is.na(padj_int) & padj_int < ALPHA)))
vs_aging <- cmp %>% group_by(intervention = panel) %>%
  summarise(aging_DEGs = nrow(aging_deg), compared = n(), spearman_rho = round(cor(aging, lfc_int, method = "spearman"), 2),
            pct_opposite = round(100 * mean(opposite), 1), DEG_in_both = sum(both_DEG),
            both_opposite = sum(both_DEG & opposite), both_same = sum(both_DEG & !opposite)) %>%
  arrange(match(intervention, c("aDR", "YMP")))
write_csv(vs_aging, file.path(OUT, "intervention_vs_aging.csv"))
vs_aging

# Scatter per intervention: grey = aging DEG only; coloured and named = DEG in aging and in the intervention
ttl <- with(vs_aging, setNames(sprintf("%s: rho = %.2f, %.0f%% opposite", intervention, spearman_rho, pct_opposite), intervention))
cmp$title <- factor(ttl[cmp$panel], ttl[c("aDR", "YMP")])
p_vs <- ggplot(cmp, aes(aging, lfc_int)) +
  geom_hline(yintercept = 0, color = "grey60") + geom_vline(xintercept = 0, color = "grey60") +
  geom_point(data = filter(cmp, !both_DEG), color = "grey70", size = 3) +
  geom_point(data = filter(cmp, both_DEG), aes(color = ifelse(opposite, "opposite to aging", "same as aging")), size = 6) +
  geom_text_repel(data = filter(cmp, both_DEG), aes(label = symbol), size = 6, seed = 42, max.overlaps = 50, fontface = "italic") +
  scale_color_manual(values = c("opposite to aging" = "#1baf7a", "same as aging" = "#a35bd6")) +
  facet_wrap(~ title) + labs(x = "log2FC aging (old vs young, shrunken)", y = "log2FC intervention (shrunken)", color = "DEG in both:") +
  guides(color = guide_legend(override.aes = list(size = 6)))
options(repr.plot.width = 20, repr.plot.height = 9.5); p_vs
save_png(p_vs, "intervention_vs_aging.png", 20, 9.5)

# Genes that are DEGs in aging and in an intervention
cmp %>% filter(both_DEG) %>% transmute(intervention = panel, symbol, aging = round(aging, 2), intervention_lfc = round(lfc_int, 2),
                                       padj_int = signif(padj_int, 2), direction = ifelse(opposite, "opposite", "same")) %>%
  arrange(intervention, direction, desc(abs(aging)))

# **Reading the comparison**
# - **aDR:** no opposition at the gene level: ρ = 0.13 and 51 % opposite sign (a coin flip would give 50 %). The 5 genes that are DEGs in both all move in the same direction as aging (*Fam163b, Etnppl, Cirbp, Rbm3, Hlf*); *Cirbp* and *Rbm3* are cold-inducible genes, possibly reflecting lower body temperature under dietary restriction.
# - **YMP:** a weak tendency towards opposition: ρ = −0.09 and 62 % opposite sign. Only 2 genes are DEGs in both, both opposite to aging (*Parp4, Hsph1*).
# - **Overall:** for most aging DEGs, the intervention's shrunken log2FC is close to 0, so neither intervention measurably changes the core aging genes (complement, microglia, MHC). The opposition seen in Section 4 is at the level of broad pathways, not of the individual aging genes.

# ====================================================================================================
# ## Summary
# **Main findings**
# - **Aging** (old vs young): 175 DEGs, mostly up; a neuroinflammation signature (complement, activated microglia, MHC class I, cathepsins, *Gfap*), with interferon and inflammatory pathways up and oxidative phosphorylation, mTORC1 and other metabolic pathways down.
# - **aDR:** 118 DEGs with small effects (glucocorticoid, circadian, immediate early genes). At the pathway level, the 6 aging sets also significant in aDR are all opposite to aging; at the gene level, aging DEGs are not opposed (ρ = 0.13, 51 % opposite).
# - **YMP:** 50 DEGs, partly meningeal genes. 9 of 11 shared aging sets opposite to aging, within a general upward shift of pathways; at the gene level a weak tendency (ρ = −0.09, 62 % opposite).
#
# **Limitations**
# - 3–5 mice per group: modest power, and fold changes of single genes in 3 vs 3 mice are uncertain.
# - The three cohorts are separate experiments with different sexes (aDR female, YMP male) and ages (intervention mice 19 months), so aDR and YMP cannot be compared with each other directly, only each with aging.
# - Bulk tissue: differences in dissection (meninges, choroid plexus in SRR23914572) can appear as expression changes.
# - A shrunken log2FC near 0 means little evidence of change, not proof of no change.
# - The comparison with aging is descriptive; no test of 'rejuvenation' is made.
#
# **Methodological decisions**
#
# | Decision | Reason |
# |---|---|
# | RSEM counts + effective lengths (tximport workflow) | corrects for gene length differences between samples |
# | Exclude SRR21355389 | haemoglobin genes 11–15× above any other aging sample (blood contamination) |
# | Keep SRR23914572 | choroid plexus genes, not blood (team decision) |
# | Pool 18 + 21 months, `~ sex + age_bin` | team decision; sex is the main axis of variation in the aging cohort |
# | Filter ≥ 10 reads in the smallest-group size; Wald test; padj < 0.05, `alpha = 0.05` | team setup |
# | apeglm-shrunken log2FC for plots and comparisons | stable fold changes for low-count genes |
# | Hallmark GSEA on the Wald statistic | uses all tested genes, not only DEGs |
# | Removed from the original notebook | correlation clustering, expression histograms, 21 vs 18 month analysis, heatmaps for aDR / YMP, GO BP, camera, permutation and binomial tests, alternative aging signatures, cell-type scores, GO of reversed genes |
# ====================================================================================================

# Software versions (reproducibility)
sessionInfo()
