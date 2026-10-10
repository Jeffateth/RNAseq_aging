# Helpers for the differential-expression notebooks (run from notebooks/).
suppressPackageStartupMessages({library(tidyverse); library(tximport); library(DESeq2); library(ggrepel)})
ROOT <- normalizePath("..")
theme_set(theme_minimal(base_size = 12))

#' The 27 runs with cohort / sex / age / group (factors ordered with the control first).
load_samples <- function() {
  read_csv(file.path(ROOT, "metadata/sra_metadata_with_group.csv"), show_col_types = FALSE) %>%
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
  dds <- if (is.matrix(x)) DESeqDataSetFromMatrix(round(x[, rownames(samples)]), samples, design)
         else DESeqDataSetFromTximport(subset_txi(x, rownames(samples)), samples, design)
  dds <- dds[rowSums(counts(dds) >= 10) >= min_n, ]
  DESeq(dds, quiet = TRUE)
}

#' Tidy results for one contrast (c(factor, treated, control)), annotated with gene symbols.
de_table <- function(dds, contrast, genes) {
  results(dds, contrast = contrast) %>% as.data.frame() %>% rownames_to_column("gene_id") %>%
    left_join(genes, by = "gene_id") %>% arrange(padj)
}

n_deg <- function(res, alpha = 0.05) res %>% filter(padj < alpha) %>% summarise(up = sum(log2FoldChange > 0), down = sum(log2FoldChange < 0))

pca_plot <- function(dds, color, shape = NULL) {
  v <- vst(dds, blind = TRUE)
  d <- plotPCA(v, intgroup = color, returnData = TRUE); pv <- round(100 * attr(d, "percentVar"))
  d[[color]] <- colData(v)[[color]]
  if (!is.null(shape)) d[[shape]] <- colData(v)[[shape]]
  ggplot(d, aes(PC1, PC2, color = .data[[color]])) +
    (if (is.null(shape)) geom_point(size = 3.5) else geom_point(aes(shape = .data[[shape]]), size = 3.5)) +
    geom_text_repel(aes(label = sub("^SRR", "", name)), size = 3, show.legend = FALSE) +
    labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"))
}

#' Genes that most distinguish one sample from the rest (vst difference) - to diagnose an outlier.
outlier_genes <- function(dds, sample, genes, n = 10) {
  v <- assay(vst(dds, blind = TRUE))
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

#' The four analyses of notebook 03 for either count source: aging 18v3 / 21v3 (sex covariate; `exclude` samples dropped),
#' aDR vs AL, YMP vs PBS. Returns a named list of tidy result tables.
run_de <- function(x, samples, genes, exclude = "SRR21355389") {
  grp <- function(g) samples %>% filter(group %in% g) %>% mutate(group = droplevels(group))
  s_age <- samples %>% filter(cohort == "aging", !rownames(.) %in% exclude) %>% mutate(group = droplevels(group))
  d_age <- make_dds(x, s_age, ~ sex + group, min_n = 3)
  d_adr <- make_dds(x, grp(c("AL", "aDR")), ~ group, min_n = 4)
  d_ymp <- make_dds(x, grp(c("PBS", "YMP")), ~ group, min_n = 3)
  list(aging_18v3 = de_table(d_age, c("group", "aging 18mo", "aging 3mo"), genes),
       aging_21v3 = de_table(d_age, c("group", "aging 21mo", "aging 3mo"), genes),
       aDR = de_table(d_adr, c("group", "aDR", "AL"), genes),
       YMP = de_table(d_ymp, c("group", "YMP", "PBS"), genes))
}

#' STAR GeneCounts matrix saved by notebook 04 (genes x runs, GENCODE vM25 versioned gene ids = kallisto's).
load_star_counts <- function(path = file.path(ROOT, "results/star/counts_unstranded.csv.gz")) {
  m <- read_csv(path, show_col_types = FALSE) %>% column_to_rownames(colnames(.)[1]) %>% as.matrix()
  storage.mode(m) <- "integer"; m
}

#' Agreement of two result tables for one contrast (a = kallisto, b = STAR, by convention).
compare_contrast <- function(a, b, alpha = 0.05) {
  j <- inner_join(a, b, by = "gene_id", suffix = c("_a", "_b")) %>% filter(!is.na(padj_a), !is.na(padj_b))
  sa <- j$padj_a < alpha; sb <- j$padj_b < alpha; sig <- sa | sb
  tibble(genes_tested = nrow(j), DEG_kallisto = sum(sa), DEG_star = sum(sb), shared = sum(sa & sb),
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
