# R helpers for notebooks 03, 04, 05 (run them from notebooks/; paths are relative to the repo root).
#
# Contents (functions grouped by purpose)
#   Load data ........ load_samples, load_txi, subset_txi, load_star_counts
#   DESeq2 ........... make_dds, de_table, n_deg, run_dds, de_tables, run_de
#   Explore .......... pca_plot, outlier_genes, volcano
#   Compare methods .. compare_de, compare_contrast, lfc_scatter          (kallisto vs STAR; outlier sensitivity)
#   Aging vs treatment overlap_with_aging, signature_score, score_test, score_plot, lfc_vs_aging_plot,
#                      msigdb_sets, run_gsea, nes_tiles, gsea_concordance, gene_table
# Terms: docs/glossary.md. Conventions: `res` = DESeq2 results table (columns gene_id, baseMean, log2FoldChange, lfcSE, stat,
# pvalue, padj, symbol); `x` = tximport list (kallisto) or gene x sample count matrix (STAR); `dds` = DESeq2 dataset.
suppressPackageStartupMessages({library(tidyverse); library(tximport); library(DESeq2); library(ggrepel)})
ROOT <- normalizePath("..")
theme_set(theme_minimal(base_size = 12))

#' The 27 runs with cohort / sex / age / group (factors ordered with the control first).
load_samples <- function() {
  read_csv(file.path(ROOT, "metadata/samples.csv"), show_col_types = FALSE) %>%
    transmute(srr = SRR, mouse_id, sex, age = age_months, group,
              cohort = if_else(BioProject == "PRJNA875066", "aging", "rejuvenation")) %>%
    mutate(group = factor(group, levels = c("aging 3mo", "aging 18mo", "aging 21mo", "AL", "aDR", "PBS", "YMP")),
           sex = factor(sex)) %>%
    column_to_rownames("srr")
}

#' kallisto -> gene-level counts (tximport). Transcript ids are GENCODE headers "tx|gene|...|symbol|..."
load_txi <- function(srr) {
  files <- setNames(file.path(ROOT, "results/kallisto", srr, "abundance.tsv"), srr)
  ids <- read_tsv(files[[1]], col_select = target_id, show_col_types = FALSE)$target_id
  f <- str_split_fixed(ids, fixed("|"), 8)
  tx2gene <- tibble(TXNAME = ids, GENEID = f[, 2])
  genes <- distinct(tibble(gene_id = f[, 2], symbol = f[, 6]))
  invisible(capture.output(txi <- suppressMessages(tximport(files, type = "kallisto", tx2gene = tx2gene))))
  list(txi = txi, genes = genes)
}

subset_txi <- function(txi, srr) {
  txi[c("abundance", "counts", "length")] <- lapply(txi[c("abundance", "counts", "length")], \(m) m[, srr, drop = FALSE])
  txi
}

#' DESeq2 dataset for a subset of samples from a tximport list (kallisto) or a gene x sample count matrix (STAR);
#' genes need >= 10 counts in >= min_n samples.
make_dds <- function(x, samples, design, min_n) {
  dds <- suppressMessages(if (is.matrix(x)) DESeqDataSetFromMatrix(round(x[, rownames(samples)]), samples, design)
                          else DESeqDataSetFromTximport(subset_txi(x, rownames(samples)), samples, design))
  dds <- dds[rowSums(counts(dds) >= 10) >= min_n, ]
  suppressMessages(DESeq(dds, quiet = TRUE))
}

#' Tidy results for one contrast (c(factor, treated, control)), annotated with gene symbols.
de_table <- function(dds, contrast, genes) {
  results(dds, contrast = contrast) %>% as.data.frame() %>% rownames_to_column("gene_id") %>%
    left_join(genes, by = "gene_id") %>% arrange(padj)
}

n_deg <- function(res, alpha = 0.05) res %>% filter(padj < alpha) %>% summarise(up = sum(log2FoldChange > 0), down = sum(log2FoldChange < 0))

pca_plot <- function(dds, color, shape = NULL) {
  v <- suppressMessages(vst(dds, blind = TRUE))
  d <- suppressMessages(plotPCA(v, intgroup = color, returnData = TRUE)); pv <- round(100 * attr(d, "percentVar"))
  d[[color]] <- colData(v)[[color]]
  if (!is.null(shape)) d[[shape]] <- colData(v)[[shape]]
  ggplot(d, aes(PC1, PC2, color = .data[[color]])) +
    (if (is.null(shape)) geom_point(size = 3.5) else geom_point(aes(shape = .data[[shape]]), size = 3.5)) +
    geom_text_repel(aes(label = sub("^SRR", "", name)), size = 3, show.legend = FALSE) +
    labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"))
}

#' Genes that most distinguish one sample from the rest (vst difference) - to diagnose an outlier.
outlier_genes <- function(dds, sample, genes, n = 10) {
  v <- assay(suppressMessages(vst(dds, blind = TRUE)))
  tibble(gene_id = rownames(v), vst_diff = v[, sample] - rowMeans(v[, colnames(v) != sample])) %>%
    left_join(genes, by = "gene_id") %>% arrange(desc(abs(vst_diff))) %>% head(n) %>% select(symbol, vst_diff)
}

#' Sensitivity of one contrast to a choice: DEG counts, overlap and LFC agreement between two result tables.
compare_de <- function(a, b, alpha = 0.05) {
  j <- inner_join(a, b, by = "gene_id", suffix = c("_a", "_b")) %>% filter(!is.na(padj_a), !is.na(padj_b))
  sig <- j$padj_a < alpha | j$padj_b < alpha
  tibble(DEG_a = sum(j$padj_a < alpha), DEG_b = sum(j$padj_b < alpha), DEG_both = sum(j$padj_a < alpha & j$padj_b < alpha),
         lfc_spearman = cor(j$log2FoldChange_a[sig], j$log2FoldChange_b[sig], method = "spearman"))
}

volcano <- function(res, title, alpha = 0.05, n_label = 12) {
  d <- res %>% filter(!is.na(padj)) %>% mutate(sig = padj < alpha, lab = if_else(row_number() <= n_label & sig, symbol, NA_character_))
  ggplot(d, aes(log2FoldChange, -log10(padj), color = sig, label = lab)) + geom_point(size = 0.8, alpha = 0.6) +
    geom_text_repel(data = filter(d, !is.na(lab)), size = 3, max.overlaps = 30, show.legend = FALSE) + scale_color_manual(values = c("grey60", "firebrick")) +
    labs(title = title, color = paste0("padj < ", alpha))
}

#' DESeq2 datasets behind the four analyses: aging cohort (`exclude` samples dropped), aDR vs AL pair, YMP vs PBS pair.
run_dds <- function(x, samples, exclude = "SRR21355389") {
  grp <- function(g) samples %>% filter(group %in% g) %>% mutate(group = droplevels(group))
  s_age <- samples %>% filter(cohort == "aging", !rownames(.) %in% exclude) %>% mutate(group = droplevels(group))
  list(aging = make_dds(x, s_age, ~ sex + group, min_n = 3),
       aDR = make_dds(x, grp(c("AL", "aDR")), ~ group, min_n = 4),
       YMP = make_dds(x, grp(c("PBS", "YMP")), ~ group, min_n = 3))
}

#' The four analyses of notebook 03 for either count source: aging 18v3 / 21v3 (sex covariate), aDR vs AL, YMP vs PBS.
#' Returns a named list of tidy result tables.
run_de <- function(x, samples, genes, exclude = "SRR21355389") de_tables(run_dds(x, samples, exclude), genes)

#' The four result tables from the datasets of run_dds().
de_tables <- function(d, genes) list(
  aging_18v3 = de_table(d$aging, c("group", "aging 18mo", "aging 3mo"), genes),
  aging_21v3 = de_table(d$aging, c("group", "aging 21mo", "aging 3mo"), genes),
  aDR = de_table(d$aDR, c("group", "aDR", "AL"), genes),
  YMP = de_table(d$YMP, c("group", "YMP", "PBS"), genes))

#' STAR GeneCounts matrix saved by notebook 01 (genes x runs, GENCODE vM25 versioned gene ids = kallisto's).
load_star_counts <- function(path = file.path(ROOT, "results/star/counts_unstranded.csv.gz")) {
  m <- read_csv(path, show_col_types = FALSE) %>% column_to_rownames(colnames(.)[1]) %>% as.matrix()
  storage.mode(m) <- "integer"; m
}

#' Agreement of two result tables for one contrast (a = kallisto, b = STAR, by convention).
compare_contrast <- function(a, b, alpha = 0.05) {
  j <- inner_join(a, b, by = "gene_id", suffix = c("_a", "_b")) %>% filter(!is.na(padj_a), !is.na(padj_b))
  sa <- j$padj_a < alpha; sb <- j$padj_b < alpha; sig <- sa | sb
  tibble(total_DEG_kallisto = sum(a$padj < alpha, na.rm = TRUE), total_DEG_star = sum(b$padj < alpha, na.rm = TRUE),   # all DEGs, as in the table above
         genes_tested_both = nrow(j), DEG_kallisto = sum(sa), DEG_star = sum(sb), shared = sum(sa & sb),   # restricted to genes tested by both
         jaccard = round(sum(sa & sb) / max(sum(sig), 1), 2),
         same_direction = round(mean(sign(j$log2FoldChange_a[sa & sb]) == sign(j$log2FoldChange_b[sa & sb])), 3),
         lfc_spearman_all = round(cor(j$log2FoldChange_a, j$log2FoldChange_b, method = "spearman"), 3),
         lfc_spearman_sig = round(cor(j$log2FoldChange_a[sig], j$log2FoldChange_b[sig], method = "spearman"), 3))
}

lfc_scatter <- function(a, b, title, alpha = 0.05) {
  j <- inner_join(a, b, by = "gene_id", suffix = c("_k", "_s")) %>% filter(!is.na(padj_k), !is.na(padj_s)) %>%
    mutate(call = case_when(padj_k < alpha & padj_s < alpha ~ "both", padj_k < alpha ~ "kallisto only",
                            padj_s < alpha ~ "STAR only", TRUE ~ "neither"),
           lab = if_else(call != "neither" & (abs(log2FoldChange_k - log2FoldChange_s) > 1), symbol_k, NA_character_))
  ggplot(j, aes(log2FoldChange_k, log2FoldChange_s, color = call)) + geom_abline(color = "grey50") +
    geom_point(data = filter(j, call == "neither"), size = 0.4, alpha = 0.3) + geom_point(data = filter(j, call != "neither"), size = 1) +
    geom_text_repel(data = filter(j, !is.na(lab)), aes(label = lab), size = 3, show.legend = FALSE, max.overlaps = 20) +
    scale_color_manual(values = c(both = "firebrick", `kallisto only` = "steelblue", `STAR only` = "darkorange", neither = "grey70")) +
    labs(title = title, x = "log2FC kallisto", y = "log2FC STAR")
}

# ---- aging vs rejuvenation ------------------------------------------------------------------
#' Overlap of an intervention's DEGs with the aging DEGs (genes tested in both): counts by direction, enrichment, log2FC agreement.
#' "reversal" = opposite sign to aging (rejuvenating); "same" = same sign (aging-like).
overlap_with_aging <- function(aging, interv, alpha = 0.05) {
  j <- inner_join(aging, interv, by = "gene_id", suffix = c("_age", "_int")) %>% filter(!is.na(padj_age), !is.na(padj_int))
  a <- j$padj_age < alpha; i <- j$padj_int < alpha; same <- sign(j$log2FoldChange_age) == sign(j$log2FoldChange_int)
  ft <- fisher.test(table(factor(a, c(FALSE, TRUE)), factor(i, c(FALSE, TRUE))), alternative = "greater")
  tibble(genes = nrow(j), aging_DEG = sum(a), intervention_DEG = sum(i), shared = sum(a & i),
         expected_by_chance = round(sum(a) * sum(i) / nrow(j), 1), shared_reversal = sum(a & i & !same), shared_same = sum(a & i & same),
         fisher_p = signif(ft$p.value, 2),
         lfc_spearman_all = round(cor(j$log2FoldChange_age, j$log2FoldChange_int, method = "spearman"), 3),
         lfc_spearman_agingDEG = round(cor(j$log2FoldChange_age[a], j$log2FoldChange_int[a], method = "spearman"), 3))
}

#' Per-sample signature score: z-score each gene across the samples of `dds`, then mean(up genes) - mean(down genes).
signature_score <- function(dds, up, down) {
  z <- t(scale(t(assay(suppressMessages(vst(dds, blind = TRUE)))))); ids <- rownames(z)
  m <- \(g) colMeans(z[intersect(g, ids), , drop = FALSE])
  tibble(srr = colnames(z), group = as.character(colData(dds)$group), score = m(up) - m(down))
}

#' Treated - control difference in score with a Welch t-test (small n: treat p-values as indicative).
score_test <- function(sc, treated, control) {
  x <- sc$score[sc$group == treated]; y <- sc$score[sc$group == control]
  tibble(treated = treated, control = control, n_t = length(x), n_c = length(y), diff = round(mean(x) - mean(y), 3),
         t_test_p = signif(t.test(x, y)$p.value, 2))
}

#' MSigDB mouse gene sets as named lists of gene symbols: Hallmark (50) and GO biological process.
msigdb_sets <- function() {
  to_list <- \(d) split(d$gene_symbol, d$gs_name)
  suppressMessages(list(hallmark = to_list(msigdbr::msigdbr(species = "Mus musculus", collection = "H")),
                        gobp = to_list(msigdbr::msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP"))))
}

#' fgsea on the DESeq2 Wald statistic (one value per gene symbol).
run_gsea <- function(res, sets, seed = 1) {
  st <- res %>% filter(!is.na(stat), !is.na(symbol), symbol != "") %>% group_by(symbol) %>% slice_max(abs(stat), n = 1, with_ties = FALSE) %>% ungroup()
  set.seed(seed)
  withCallingHandlers(fgsea::fgsea(sets, setNames(st$stat, st$symbol), minSize = 15, maxSize = 500),
                      warning = \(w) if (grepl("ties in the preranked", conditionMessage(w))) invokeRestart("muffleWarning")) %>%   # ~0.1 % tied stats: harmless
    as_tibble() %>% select(-leadingEdge)
}

#' Intervention log2FC against aging log2FC (genes tested in both); negative slope among aging DEGs = reversal.
lfc_vs_aging_plot <- function(aging, interv, title, alpha = 0.05) {
  j <- inner_join(aging, interv, by = "gene_id", suffix = c("_age", "_int")) %>% filter(!is.na(padj_age), !is.na(padj_int)) %>%
    mutate(call = case_when(padj_age < alpha & padj_int < alpha ~ "both", padj_age < alpha ~ "aging only",
                            padj_int < alpha ~ "intervention only", TRUE ~ "neither"))
  a <- j$padj_age < alpha
  ggplot(j, aes(log2FoldChange_age, log2FoldChange_int, color = call)) +
    geom_hline(yintercept = 0, color = "grey85") + geom_vline(xintercept = 0, color = "grey85") +
    geom_point(data = filter(j, call == "neither"), size = 0.4, alpha = 0.3) + geom_point(data = filter(j, call != "neither"), size = 1.3) +
    scale_color_manual(values = c(both = "firebrick", `aging only` = "steelblue", `intervention only` = "darkorange", neither = "grey70")) +
    labs(title = title, subtitle = sprintf("Spearman among aging DEGs = %.2f", cor(j$log2FoldChange_age[a], j$log2FoldChange_int[a], method = "spearman")),
         x = "log2FC aging (21 vs 3 mo)", y = "log2FC intervention")
}

#' Per-sample signature scores for the three datasets (list with aging / aDR / YMP data frames from signature_score()).
score_plot <- function(sc) {
  d <- bind_rows(aging = sc$aging, aDR = sc$aDR, YMP = sc$YMP, .id = "panel") %>%
    mutate(group = factor(group, c("aging 3mo", "aging 18mo", "aging 21mo", "AL", "aDR", "PBS", "YMP")))
  ggplot(d, aes(group, score, color = group)) + geom_boxplot(outlier.shape = NA, fill = NA) + geom_jitter(width = 0.1, size = 2) +
    facet_wrap(~ panel, scales = "free_x") + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1)) + labs(x = NULL, y = "aging-signature score")
}

#' Hallmark (or other) NES tiles for aging vs the interventions; * = padj < 0.05; only sets significant in at least one contrast.
nes_tiles <- function(g, alpha = 0.05) {
  d <- bind_rows(g, .id = "contrast") %>% mutate(pathway = sub("^HALLMARK_", "", pathway))
  keep <- unique(d$pathway[d$padj < alpha]); ord <- d %>% filter(contrast == "aging_21v3") %>% arrange(NES) %>% pull(pathway)
  d <- d %>% filter(pathway %in% keep) %>% mutate(pathway = factor(pathway, ord), contrast = factor(contrast, c("aging_21v3", "aDR", "YMP")))
  ggplot(d, aes(contrast, pathway, fill = NES)) + geom_tile(color = "white") + geom_text(aes(label = ifelse(padj < alpha, "*", "")), size = 5) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") + labs(x = NULL, y = NULL)
}

#' Pathways significant in aging: how many does each intervention reverse (opposite NES sign) or mimic (same sign), both padj < alpha.
gsea_concordance <- function(g, alpha = 0.05) {
  bind_rows(lapply(c("aDR", "YMP"), \(n) {
    j <- inner_join(g$aging_21v3 %>% select(pathway, NES_age = NES, padj_age = padj), g[[n]] %>% select(pathway, NES, padj), by = "pathway")
    s <- filter(j, padj_age < alpha)
    tibble(intervention = n, aging_sets = nrow(s), significant_in_both = sum(s$padj < alpha),
           reversed = sum(s$padj < alpha & sign(s$NES) != sign(s$NES_age)), same_direction = sum(s$padj < alpha & sign(s$NES) == sign(s$NES_age)),
           NES_spearman_all_sets = round(cor(j$NES_age, j$NES, method = "spearman"), 2))
  }))
}

#' log2FC / padj of named genes in the aging and intervention contrasts.
gene_table <- function(res, symbols) {
  pick <- \(r, tag) r %>% filter(symbol %in% symbols) %>% transmute(symbol, !!paste0("lfc_", tag) := round(log2FoldChange, 2), !!paste0("padj_", tag) := signif(padj, 2))
  Reduce(\(a, b) left_join(a, b, by = "symbol"), list(pick(res$aging_21v3, "aging"), pick(res$aDR, "aDR"), pick(res$YMP, "YMP")))
}
