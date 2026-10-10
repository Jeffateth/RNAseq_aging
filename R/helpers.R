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

#' DESeq2 dataset for a subset of samples; genes need >= 10 counts in >= min_n samples.
make_dds <- function(txi, samples, design, min_n) {
  dds <- DESeqDataSetFromTximport(subset_txi(txi, rownames(samples)), samples, design)
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
