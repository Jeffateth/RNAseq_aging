# Lea's analysis on the team setup (proposed baseline). Exported from teammates/lea/Script.ipynb; text cells as comments.

# ====================================================================================================
# # Aging and rejuvenating interventions in the mouse cerebellum
#
# > **Proposed common baseline for review, not final.** Lea's notebook reduced to the team's agreed setup. Results, changes from Lea's original, open decisions and slide check: `README.md`.
#
# **Question:** how does the cerebellum change with age, what do dietary restriction (aDR) and young plasma (YMP) change, and does either move the aging changes back?
#
# - **Data:** RSEM counts (GENCODE M39) of 27 mice from Hahn et al. (2023), in `../../for_teammates/`.
# - **Run:** `pixi run jupyter nbconvert --to notebook --execute --inplace Script.ipynb` (outputs in `output/`; the first run downloads MSigDB 2026.1 once).
# - **Comparisons**, each within its own cohort: aging (old vs young), aDR vs AL, YMP vs PBS.
#
# | Section | Question |
# |---|---|
# | 1 · Import and QC | Which samples go in? |
# | 2 · PCA | What drives the main variation? |
# | 3 · Differential expression | Which genes change? |
# | 4 · Gene sets (GSEA) | Which processes shift? |
# | 5 · Aging vs interventions | Do aDR or YMP oppose aging? |
# ====================================================================================================

# Packages, figure style and output folder
suppressPackageStartupMessages({
  library(DESeq2); library(apeglm); library(fgsea); library(msigdbr); library(ggrepel)
  library(tidyverse)   # loaded last, so dplyr's count/select/filter are not masked
})
options(jupyter.rich_display = FALSE, width = 140)
# Large text and tight margins in every figure
theme_set(theme_minimal(base_size = 40) +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "#333333"), plot.title = element_text(size = rel(1)),
        strip.text = element_text(size = rel(1)), legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0.2, "lines"), legend.key.size = unit(2, "lines"), legend.spacing = unit(0.2, "lines"),
        panel.spacing = unit(1.5, "lines"), plot.margin = margin(2, 6, 2, 2), axis.title = element_text(size = rel(1))))
LAB <- 12   # text labels inside plots (mm, about 34 pt)
DATA <- "../../for_teammates"   # shared RSEM tables (in git)
OUT  <- "output"; dir.create(OUT, showWarnings = FALSE)
ALPHA <- 0.05                   # DEG and gene set cutoff (padj), also DESeq2's results(alpha = )
COL <- c(aging = "#2a78d6", aDR = "#eb6834", YMP = "#1baf7a", up = "#e34948", down = "#2a78d6", ns = "#cfcfcf")
save_png <- function(p, name, w, h) ggsave(file.path(OUT, name), p, width = w, height = h, dpi = 150, bg = "white")

# ====================================================================================================
# ## 1 · Import and quality control
# RSEM expected counts, TPM and effective gene lengths. DESeq2 gets counts plus lengths (`tximport` workflow), so it corrects for sequencing depth and per-sample gene length. The notebook stops if any length is missing.
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
# - **Excluded: SRR21355389** (21 months): haemoglobin 11–15× above any other aging sample (table: `x_vs_max`), i.e. blood in the tissue, which would appear as false aging changes.
# - **Kept: SRR23914572** (AL): its unusual genes are choroid plexus, not blood (team decision).
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
# ## 2 · PCA
# Variance-stabilised counts (`vst`) of genes with ≥ 10 reads in ≥ 3 samples, 500 most variable genes. Each point is a mouse; close = similar profiles.
# ====================================================================================================

# VST of all 26 samples, then PCA (top 500 genes); colour and shape = cohort
vst_of <- function(srr, min_n = 3) {
  x <- txi; x[c("abundance", "counts", "length")] <- lapply(x[c("abundance", "counts", "length")], \(m) m[, srr])
  d <- suppressMessages(DESeqDataSetFromTximport(x, meta[srr, ], ~ 1))
  suppressMessages(vst(d[rowSums(counts(d) >= 10) >= min_n, ], blind = TRUE))
}
pca_df <- function(v) {
  p <- suppressMessages(plotPCA(v, intgroup = "group", returnData = TRUE, ntop = 500))
  list(d = cbind(p, meta[rownames(p), c("experiment", "sex")]), pv = round(100 * attr(p, "percentVar")))
}
v_all <- vst_of(rownames(meta))
p <- pca_df(v_all)
p$d$experiment <- factor(p$d$experiment, c("aging", "aDR", "YMP"))
p_all <- ggplot(p$d, aes(PC1, PC2, color = experiment, shape = experiment)) + geom_point(size = 10, alpha = 0.85) +
  scale_color_manual(values = COL) + labs(x = paste0("PC1 (", p$pv[1], "%)"), y = paste0("PC2 (", p$pv[2], "%)"),
                                          color = NULL, shape = NULL) + guides(color = guide_legend(override.aes = list(size = 10)))
options(repr.plot.width = 10, repr.plot.height = 8); p_all
save_png(p_all, "pca_all.png", 10, 8)

# **Reading:** PC1 (41 %) separates the aging cohort from the intervention cohorts, PC2 (13 %) aDR from YMP. Experiments differ more than conditions, so each comparison stays within its cohort.

# ====================================================================================================
# ## 3 · Differential expression (DESeq2)
# | Choice | Setting (agreed) | Why |
# |---|---|---|
# | Aging | old (18 + 21 mo, n = 7) vs young (3 mo, n = 4), `~ sex + age_bin` | both sexes in this cohort: the sex term keeps sex differences out of the age test |
# | aDR, YMP | aDR (5) vs AL (4), YMP (3) vs PBS (3), `~ condition` | single-sex cohorts |
# | Gene filter | ≥ 10 reads in ≥ smallest-group size (4 / 4 / 3) | a gene can be expressed in one whole group |
# | Test, DEG | Wald, `alpha = 0.05`; DEG = padj < 0.05, no fold-change cutoff | |
# | Fold changes in plots | apeglm-shrunken | tames noisy low-count genes; p-values unchanged |
#
# - **Pooling 18 + 21 months:** in Lea's original analysis only 16 genes differed between 21 and 18 months (vs 143–149 against 3 months); the heatmap shows them alike.
# - **Tested genes:** genes passing the filter (they get a Wald statistic). Some get padj = NA (too few reads, DESeq2's independent filtering): not DEGs, but still ranked in GSEA.
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

meta$age_bin   <- ifelse(meta$experiment != "aging", NA, ifelse(meta$age_months == 3, "young", "old"))   # aging cohort only
meta$condition <- as.character(meta$group)
srr <- split(rownames(meta), meta$experiment)
dds <- list(
  aging = fit_de(srr$aging, ~ sex + age_bin, "age_bin",   c("young", "old"), min_n = 4),
  aDR   = fit_de(srr$aDR,   ~ condition,     "condition", c("AL", "aDR"),    min_n = 4),
  YMP   = fit_de(srr$YMP,   ~ condition,     "condition", c("PBS", "YMP"),   min_n = 3))
res <- lapply(dds, de_table)
sapply(dds, \(d) tail(resultsNames(d), 1))   # tested coefficient per comparison

# DEG summary. reference = counts from the team's reference run (results/common), to check that the setup is reproduced
deg_summary <- imap_dfr(res, \(r, n) tibble(comparison = n, tested = sum(!is.na(r$stat)), with_padj = sum(!is.na(r$padj)),
  DEGs = sum(r$padj < ALPHA, na.rm = TRUE), up = sum(r$padj < ALPHA & r$log2FoldChange > 0, na.rm = TRUE),
  down = sum(r$padj < ALPHA & r$log2FoldChange < 0, na.rm = TRUE))) %>%
  mutate(reference = c(175, 118, 50))
deg_summary

# Volcano plots: x = shrunken log2FC, colour = DEG (padj < 0.05) up / down; the 6 most significant DEGs are named.
# x axis limited to -5..5: the few genes beyond are drawn at the edge (listed below the plot)
vd <- bind_rows(res, .id = "comparison") %>% filter(!is.na(padj)) %>%   # genes with padj only (y axis)
  mutate(comparison = factor(comparison, names(res)),
         call = case_when(padj < ALPHA & log2FoldChange > 0 ~ "up", padj < ALPHA ~ "down", TRUE ~ "ns"),
         title = paste0(comparison, ": ", ave(call == "up", comparison, FUN = sum), " up, ", ave(call == "down", comparison, FUN = sum), " down"),
         title = fct_reorder(title, as.integer(comparison)))
top <- vd %>% filter(call != "ns") %>% group_by(comparison) %>% slice_min(padj, n = 6, with_ties = FALSE)
p_volc <- ggplot(vd, aes(lfc_shrunk, -log10(padj), color = call)) + geom_point(size = 2.5) +
  geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", color = "grey50") +
  geom_text_repel(data = top, aes(label = symbol), color = "black", size = LAB, fontface = "italic", seed = 42, max.overlaps = 30,
                  box.padding = 0.6, min.segment.length = 0) +
  scale_color_manual(values = COL, breaks = c("up", "down"), labels = c("up", "down")) +
  scale_x_continuous(limits = c(-5, 5), oob = scales::squish) +
  facet_wrap(~ title, scales = "free_y") + labs(x = "log2 fold change (shrunken)", y = "-log10 padj", color = NULL) +
  guides(color = guide_legend(override.aes = list(size = 9)))
options(repr.plot.width = 26, repr.plot.height = 10); p_volc
save_png(p_volc, "volcano.png", 26, 10)
vd %>% filter(call != "ns", abs(lfc_shrunk) > 5) %>% transmute(comparison, symbol, lfc_shrunk = round(lfc_shrunk, 1), padj = signif(padj, 2))

# Top 10 DEGs per comparison (shrunken log2FC; positive = higher in old / treated mice)
bind_rows(res, .id = "comparison") %>% filter(padj < ALPHA) %>% group_by(comparison) %>%
  slice_min(padj, n = 10, with_ties = FALSE) %>%
  summarise(top10 = paste0(symbol, " (", sprintf("%+.1f", lfc_shrunk), ")", collapse = ", ")) %>%
  mutate(comparison = factor(comparison, names(res))) %>% arrange(comparison)

# Heatmap of the 30 most significant aging DEGs (one row per gene, one column per mouse): VST of the aging samples, z-score per gene; columns = mice by age and sex
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
  scale_y_discrete(expand = c(0, 0)) + labs(x = NULL, y = NULL, fill = "z") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(0.8)), axis.text.y = element_text(face = "italic", size = rel(0.7)),
        panel.grid = element_blank(), legend.position = "right", legend.key.height = unit(4, "lines"))
options(repr.plot.width = 12, repr.plot.height = 16); p_heat
save_png(p_heat, "heatmap_aging.png", 12, 16)

# Save the three DE tables (all tested genes; padj NA = removed by independent filtering or a count outlier)
iwalk(res, \(r, n) write_csv(r, file.path(OUT, paste0("de_", n, ".csv"))))
write_csv(deg_summary, file.path(OUT, "deg_summary.csv"))

# **Reading**
# - DEG counts match the team reference: aging 175, aDR 118, YMP 50.
# - **Volcano:** right = higher in old / treated; dashed line = padj 0.05; genes beyond ±5 drawn at the edge (listed).
# - **Aging (138 up, 37 down)**, by category (table from Lea's original):
#
#   | Category | Genes |
#   |---|---|
#   | Activated microglia | *Clec7a, Cst7, Trem2, Tyrobp, Lpl, Gpnmb, Lgals3, Cd68, Lilrb4a* |
#   | Complement | *C4b, C1qa, C1qb, C1qc* |
#   | MHC class I | *H2-K1, H2-Q6, B2m* |
#   | Lysosomal proteases | *Ctss, Ctsd, Ctsz* |
#   | Reactive astrocytes | *Gfap* |
# - **Heatmap** (z-score per gene across the 11 aging mice): young mice are low for all 30 genes; 18 and 21 months look alike.
# - **aDR (46 up, 72 down):** small changes; glucocorticoid / circadian (*Tsc22d3, Hlf, Ciart*) and immediate early genes (*Nr4a2, Nr4a3*). *Cdr1* is on/off (0 reads in all AL mice): treat with caution.
# - **YMP (22 up, 28 down):** large fold changes rest on 3 vs 3 mice; top genes *Ptgds, Ogn, Slc47a1* are meningeal, so dissection may contribute.

# ====================================================================================================
# ## 4 · Gene sets (GSEA)
# - **Method:** all tested genes ranked by the Wald statistic; `fgsea` (sets of 15–500 genes, seed 42). **NES** > 0 = set higher in old / treated mice; significant = padj < 0.05.
# - **Hallmark** (50 curated sets) is the main result. **GO BP** (~3,700 terms) is a supplementary table, used only to name a process behind a theme already seen (e.g. myelin, circadian), since thousands of overlapping terms will always contain chance hits.
# ====================================================================================================

# Hallmark GSEA for each comparison: all 50 sets tested; the table counts the significant ones (padj < 0.05)
hallmark <- msigdbr(species = "Mus musculus", db_species = "MM", collection = "MH")
hallmark <- split(hallmark$gene_symbol, hallmark$gs_name)
cat("msigdbr", as.character(packageVersion("msigdbr")), "|", length(hallmark), "Hallmark sets\n")
run_gsea <- function(r, sets = hallmark) {
  r <- r %>% filter(!is.na(stat), !is.na(symbol), symbol != "") %>% group_by(symbol) %>%
    slice_max(abs(stat), n = 1, with_ties = FALSE) %>% ungroup()
  set.seed(42)
  suppressWarnings(fgsea(sets, setNames(r$stat, r$symbol), minSize = 15, maxSize = 500)) %>%
    as_tibble() %>% select(pathway, NES, pval, padj, size) %>% arrange(padj)
}
gsea <- bind_rows(lapply(res, run_gsea), .id = "comparison")
write_csv(gsea, file.path(OUT, "gsea_hallmark.csv"))
gsea %>% group_by(comparison = factor(comparison, names(res))) %>%
  summarise(significant = sum(padj < ALPHA), up = sum(padj < ALPHA & NES > 0), down = sum(padj < ALPHA & NES < 0))

# Supplement: GO BP GSEA, same ranking (Wald statistic), seed 42 and set sizes as Hallmark; all terms saved
gobp <- msigdbr(species = "Mus musculus", db_species = "MM", collection = "M5", subcollection = "GO:BP")
gobp <- split(gobp$gene_symbol, gobp$gs_name)
gsea_go <- bind_rows(lapply(res, run_gsea, sets = gobp), .id = "comparison")
write_csv(gsea_go, file.path(OUT, "gsea_gobp.csv"))
gsea_go %>% group_by(comparison = factor(comparison, names(res))) %>%
  summarise(tested = n(), significant = sum(padj < ALPHA), up = sum(padj < ALPHA & NES > 0), down = sum(padj < ALPHA & NES < 0))

# GO BP terms used on the current slides: NES per comparison (* = padj < 0.05). One place to take slide values from
slide_terms <- c("GOBP_COMPLEMENT_ACTIVATION", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION", "GOBP_MACROPHAGE_ACTIVATION",
                 "GOBP_NEUROINFLAMMATORY_RESPONSE", "GOBP_SYNAPSE_PRUNING", "GOBP_ASTROCYTE_ACTIVATION",
                 "GOBP_OXIDATIVE_PHOSPHORYLATION", "GOBP_STEROL_BIOSYNTHETIC_PROCESS", "GOBP_CIRCADIAN_RHYTHM",
                 "GOBP_RESPONSE_TO_STARVATION", "GOBP_SPROUTING_ANGIOGENESIS", "GOBP_RESPONSE_TO_TYPE_I_INTERFERON",
                 "GOBP_ENSHEATHMENT_OF_NEURONS", "GOBP_MYELIN_ASSEMBLY", "GOBP_OLIGODENDROCYTE_DIFFERENTIATION")
gsea_go %>% filter(pathway %in% slide_terms) %>%
  mutate(value = sprintf("%+.2f%s", NES, ifelse(padj < ALPHA, "*", "")),
         term = tolower(gsub("_", " ", sub("^GOBP_", "", pathway)))) %>%
  select(term, comparison, value) %>% pivot_wider(names_from = comparison, values_from = value) %>%
  select(term, aging, aDR, YMP) %>% arrange(match(term, tolower(gsub("_", " ", sub("^GOBP_", "", slide_terms)))))

# One figure: the Hallmark sets that the aging GSEA found significant (padj < 0.05; 15 of 50), ordered by aging NES,
# with the NES of the same sets in aging, aDR and YMP.
# Colour = NES (red = higher in old / treated mice, blue = lower); solid tile with * = padj < 0.05; faded tile = not significant
nice <- \(p) str_to_sentence(gsub("_", " ", sub("^HALLMARK_", "", p))) %>%
  str_replace_all(c("^Tnfa signaling via nfkb" = "TNFa signaling via NF-kB", "^Il6 jak stat3" = "IL6-JAK-STAT3", "^Il2 stat5" = "IL2-STAT5",
                    "^Mtorc1" = "mTORC1", "^Myc" = "MYC", "^E2f" = "E2F", "^G2m" = "G2M", "^Kras" = "KRAS", "^P53" = "p53", "^Uv" = "UV",
                    "^Tgf beta" = "TGF-beta", "^Pi3k akt mtor" = "PI3K-AKT-mTOR", "^Dna" = "DNA", "^Wnt beta catenin" = "WNT-beta-catenin"))
ord <- gsea %>% filter(comparison == "aging", padj < ALPHA) %>% arrange(NES) %>% pull(pathway)
gd <- gsea %>% filter(pathway %in% ord) %>%
  mutate(set = factor(nice(pathway), nice(ord)), comparison = factor(comparison, names(res)),
         sig = factor(ifelse(padj < ALPHA, "padj < 0.05", "not significant"), c("padj < 0.05", "not significant")))
lim <- max(abs(gd$NES))
p_gsea <- ggplot(gd, aes(comparison, set, fill = NES, alpha = sig)) + geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = ifelse(padj < ALPHA, "*", "")), size = 18, vjust = 0.75, alpha = 1) +
  scale_fill_gradient2(low = COL[["down"]], mid = "#f4f4f4", high = COL[["up"]], limits = c(-lim, lim)) +
  scale_alpha_manual(values = c("padj < 0.05" = 1, "not significant" = 0.3), name = NULL) +
  scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + labs(x = NULL, y = NULL) +
  theme(panel.grid = element_blank(), legend.position = "right", legend.key.height = unit(3, "lines")) +
  guides(alpha = guide_legend(override.aes = list(fill = COL[["up"]])), fill = guide_colorbar(order = 1))
options(repr.plot.width = 17, repr.plot.height = 12); p_gsea
save_png(p_gsea, "gsea_hallmark.png", 17, 12)

# Descriptive count (Lea's 5b step 3 definition: reversed = significant in both, opposite NES sign): of the 15 aging sets, how many are also significant in each intervention, in the opposite or the same direction
aging_sets <- gsea %>% filter(comparison == "aging", padj < ALPHA) %>% select(pathway, NES_aging = NES)
gsea %>% filter(comparison != "aging", padj < ALPHA) %>% inner_join(aging_sets, by = "pathway") %>%
  group_by(intervention = comparison) %>%
  summarise(aging_sets = nrow(aging_sets), significant_in_both = n(),
            opposite = sum(sign(NES) != sign(NES_aging)), same = sum(sign(NES) == sign(NES_aging)),
            opposite_sets = paste(nice(pathway[sign(NES) != sign(NES_aging)]), collapse = ", "),
            same_sets = paste(nice(pathway[sign(NES) == sign(NES_aging)]), collapse = ", ")) %>%
  arrange(match(intervention, c("aDR", "YMP")))

# **Reading**
# - **Hallmark:** 15 of 50 sets significant in aging (up: interferon, complement, inflammation, allograft rejection, coagulation, KRAS; down: OXPHOS, mTORC1, fatty acid, MYC, TNFα, glycolysis, hypoxia, androgen), 16 in aDR, 24 in YMP (22 up). The figure shows the 15 aging sets; faded = not significant in that comparison.
# - **Reversed set** (Lea's definition): significant in both, opposite NES. The intervention fold change is treated vs untreated *old* mice, i.e. the intervention's own effect, so reversal means an opposite sign. "Still up but less" cannot be tested: it needs young mice in the intervention cohorts.
# - **aDR:** all 6 shared aging sets opposite (complement, coagulation, KRAS down; OXPHOS, MYC, TNFα up). **YMP:** 9 of 11 opposite, 2 same (complement, coagulation); since almost all YMP sets go up, this is partly a general shift.
# - **GO BP terms on the slides:** see the table above (* = significant). Myelin terms go down with aDR and up with YMP but do not change with age, so they are not a reversal; they may also reflect white-matter content of the punch.
# - Counts are descriptive: gene sets share genes.

# ====================================================================================================
# ## 5 · Aging vs interventions
# For the 175 aging DEGs: aging vs intervention shrunken log2FC, one panel per intervention (upper-left / lower-right = opposite to aging). Described by Spearman ρ, % opposite sign and median |log2FC|; genes that are also intervention DEGs are highlighted.
# - **Only fold changes**, because the cohorts are separate experiments.
# - **No significance test:** aging genes are correlated, so a binomial or correlation p-value would overstate the evidence. Near-zero fold changes have near-random signs, so read % opposite together with the median |log2FC|.
# ====================================================================================================

# Compare each intervention with aging over the aging DEGs
aging_deg <- res$aging %>% filter(padj < ALPHA) %>% select(gene_id, symbol, aging = lfc_shrunk)
cmp <- bind_rows(lapply(c("aDR", "YMP"), \(n) aging_deg %>%
  inner_join(res[[n]] %>% select(gene_id, lfc_int = lfc_shrunk, padj_int = padj), by = "gene_id") %>%
  filter(!is.na(lfc_int)) %>%
  mutate(panel = n, opposite = sign(aging) != sign(lfc_int), both_DEG = !is.na(padj_int) & padj_int < ALPHA)))
vs_aging <- cmp %>% group_by(intervention = panel) %>%
  summarise(aging_DEGs = nrow(aging_deg), compared = n(), spearman_rho = round(cor(aging, lfc_int, method = "spearman"), 2),
            pct_opposite = round(100 * mean(opposite), 1),
            median_abs_lfc_aging = round(median(abs(aging)), 2), median_abs_lfc_int = round(median(abs(lfc_int)), 2), DEG_in_both = sum(both_DEG),
            both_opposite = sum(both_DEG & opposite), both_same = sum(both_DEG & !opposite)) %>%
  arrange(match(intervention, c("aDR", "YMP")))
write_csv(vs_aging, file.path(OUT, "intervention_vs_aging.csv"))
vs_aging

# Scatter per intervention: grey = aging DEG only; coloured and named = DEG in aging and in the intervention
ttl <- with(vs_aging, setNames(sprintf("%s: rho = %.2f, %.0f%% opposite", intervention, spearman_rho, pct_opposite), intervention))
cmp$title <- factor(ttl[cmp$panel], ttl[c("aDR", "YMP")])
p_vs <- ggplot(cmp, aes(aging, lfc_int)) +
  geom_hline(yintercept = 0, color = "grey60") + geom_vline(xintercept = 0, color = "grey60") +
  geom_point(data = filter(cmp, !both_DEG), color = "grey70", size = 4) +
  geom_point(data = filter(cmp, both_DEG), aes(color = ifelse(opposite, "opposite to aging", "same as aging")), size = 9) +
  geom_text_repel(data = filter(cmp, both_DEG), aes(label = symbol), size = LAB, seed = 42, max.overlaps = 50, fontface = "italic",
                  box.padding = 0.8, min.segment.length = 0) +
  scale_color_manual(values = c("opposite to aging" = "#1baf7a", "same as aging" = "#a35bd6")) +
  facet_wrap(~ title) + labs(x = "log2FC aging (shrunken)", y = "log2FC intervention", color = "DEG in both:") +
  guides(color = guide_legend(override.aes = list(size = 9)))
options(repr.plot.width = 22, repr.plot.height = 10); p_vs
save_png(p_vs, "intervention_vs_aging.png", 22, 10)

# Genes that are DEGs in aging and in an intervention
cmp %>% filter(both_DEG) %>% transmute(intervention = panel, symbol, aging = round(aging, 2), intervention_lfc = round(lfc_int, 2),
                                       padj_int = signif(padj_int, 2), direction = ifelse(opposite, "opposite", "same")) %>%
  arrange(intervention, direction, desc(abs(aging)))

# **Reading**
# - Aging DEGs change by a median |log2FC| of 0.70 with age but **0.01** after aDR or YMP: neither intervention measurably moves them (flat band at y = 0).
# - **aDR:** ρ = 0.13, 51 % opposite (no opposition); the 5 genes significant in both go the *same* way as aging (*Fam163b, Etnppl, Cirbp, Rbm3, Hlf*; *Cirbp*, *Rbm3* are cold-inducible).
# - **YMP:** ρ = −0.09, 62 % opposite, mostly noise at this effect size; 2 genes significant in both, both opposite (*Parp4, Hsph1*).
# - **Pathway vs gene level:** both interventions oppose some aging pathways (Section 4), but not the aging genes themselves.
# - **Confounding:** the aging signature comes from another experiment with both sexes; aDR is female-only and YMP male-only, so batch and sex-specific aging can distort the comparison. Dissection differences can mimic changes in all comparisons.

# ====================================================================================================
# ## Summary
# - **Aging:** 175 DEGs; neuroinflammation up (complement, microglia, MHC I, *Gfap*), metabolism down (OXPHOS, mTORC1).
# - **aDR:** 118 DEGs, small effects; opposes 6 aging pathways, not the aging genes.
# - **YMP:** 50 DEGs, partly meningeal; opposes 9 of 11 shared aging pathways within a general upward shift, not the aging genes.
# - **Limitations:** 3–5 mice per group; separate experiments with different sexes and ages; bulk tissue; a near-zero fold change is not proof of no change; the comparison with aging is descriptive.
# ====================================================================================================

# Software versions (reproducibility)
sessionInfo()
