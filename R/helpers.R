# R helpers for notebooks 03, 04, 05 (run them from notebooks/; paths are relative to the repo root).
#
# Contents (functions grouped by purpose)
#   Load data ........ load_samples, load_genes, load_rsem, load_star_counts, load_txi (kallisto), subset_txi
#   DESeq2 ........... make_dds, de_table, n_deg, run_dds, de_tables
#   Explore .......... pca_plot, sample_tree, outlier_genes, deg_patterns, deg_pattern_plot
#   Compare methods .. compare_de, compare_contrast, lfc_scatter          (any two methods; outlier sensitivity)
#                      counts_per_sample, method_agreement, single_method_calls, method_evidence (notebook 04)
#                      tpm_variation, tpm_variation_table, tpm_variation_plot (technical vs biological variation)
#   Aging vs treatment overlap_with_aging, signature_score, score_test, score_plot, lfc_vs_aging_plot,
#                      msigdb_sets, run_gsea, nes_tiles, gsea_concordance, gene_table
# Terms: docs/glossary.md. Conventions: `res` = DESeq2 results table (columns gene_id, baseMean, log2FoldChange, lfcSE, stat,
# pvalue, padj, symbol); `x` = tximport list (kallisto, RSEM) or gene x sample count matrix (STAR); `dds` = DESeq2 dataset.
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

#' Gene id -> symbol table of GENCODE M39 (metadata/genes_M39.csv, built by scripts/build_genes.py).
load_genes <- function() read_csv(file.path(ROOT, "metadata/genes_M39.csv"), show_col_types = FALSE)

#' kallisto -> gene-level counts (tximport). Transcript ids are GENCODE headers "tx|gene|...|symbol|..."
load_txi <- function(srr, dir = file.path(ROOT, "results/kallisto")) {
  files <- setNames(file.path(dir, srr, "abundance.tsv"), srr)
  ids <- read_tsv(files[[1]], col_select = target_id, show_col_types = FALSE)$target_id
  f <- str_split_fixed(ids, fixed("|"), 8)
  tx2gene <- tibble(TXNAME = ids, GENEID = f[, 2])
  genes <- distinct(tibble(gene_id = f[, 2], symbol = f[, 6]))
  invisible(capture.output(txi <- suppressMessages(tximport(files, type = "kallisto", tx2gene = tx2gene, dropInfReps = TRUE))))
  list(txi = txi, genes = genes)
}

#' RSEM -> gene-level counts (tximport reads RSEM's gene results directly: expected counts, TPM, effective length).
#' Genes whose effective length is 0 (not expressed) get length 1, as tximport requires lengths > 0.
load_rsem <- function(srr, dir = file.path(ROOT, "results/rsem")) {
  files <- setNames(file.path(dir, srr, paste0(srr, ".genes.results")), srr)
  invisible(capture.output(txi <- suppressMessages(tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE))))
  txi$length[txi$length <= 0] <- 1
  txi
}

subset_txi <- function(txi, srr) {
  txi[c("abundance", "counts", "length")] <- lapply(txi[c("abundance", "counts", "length")], \(m) m[, srr, drop = FALSE])
  txi
}

#' DESeq2 dataset for a subset of samples from a tximport list (kallisto, RSEM) or a gene x sample count matrix (STAR);
#' genes need >= 10 counts in >= min_n samples.
make_dds <- function(x, samples, design, min_n) {
  dds <- suppressMessages(if (is.matrix(x)) DESeqDataSetFromMatrix(round(x[, rownames(samples)]), samples, design)
                          else DESeqDataSetFromTximport(subset_txi(x, rownames(samples)), samples, design))
  dds <- dds[rowSums(counts(dds) >= 10) >= min_n, ]
  withCallingHandlers(suppressMessages(DESeq(dds, quiet = TRUE)),
                      warning = \(w) if (grepl("design is ~ 1", conditionMessage(w))) invokeRestart("muffleWarning"))   # intercept-only model of the PCA is intended
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

#' Hierarchical clustering of the samples (tutorial 3-3): distance 1 - Pearson r of the vst values (blind to groups), average linkage.
#' Labels: group, sex (f/m), run number.
sample_tree <- function(dds) {
  v <- assay(suppressMessages(vst(dds, blind = TRUE)))
  hc <- hclust(as.dist(1 - cor(v)), method = "average")
  cd <- colData(dds)
  plot(hc, labels = paste(cd$group, substr(cd$sex, 1, 1), sub("^SRR", "", colnames(v))), hang = -1, main = NULL, xlab = "", sub = "", ylab = "1 - Pearson r")
  invisible(hc)
}

#' Aging DEGs (padj < alpha in 18 vs 3 or 21 vs 3 months) grouped by their course over age (tutorial 3-5). Profile per gene = log2FC vs
#' 3 months at 3, 18, 21 months (0, then DESeq2's estimates, so adjusted for sex); distance 1 - Pearson r of the profiles (shape, not size),
#' average linkage, cut into k clusters.
deg_patterns <- function(age, k = 4, alpha = 0.05) {
  ids <- union(filter(age$`18v3`, padj < alpha)$gene_id, filter(age$`21v3`, padj < alpha)$gene_id)
  lfc <- \(t) unname(setNames(t$log2FoldChange, t$gene_id)[ids])
  m <- cbind(0, lfc(age$`18v3`), lfc(age$`21v3`))
  tibble(gene_id = ids, cluster = unname(cutree(hclust(as.dist(1 - cor(t(m))), method = "average"), k)), lfc_18v3 = m[, 2], lfc_21v3 = m[, 3]) %>%
    left_join(select(age$`18v3`, gene_id, symbol), by = "gene_id")
}

#' One panel per cluster: each gene's log2FC vs 3 months over age (grey) and the cluster median (blue).
deg_pattern_plot <- function(p) {
  d <- p %>% mutate(m3 = 0) %>% pivot_longer(c(m3, lfc_18v3, lfc_21v3), values_to = "lfc") %>%
    mutate(months = c(m3 = 3, lfc_18v3 = 18, lfc_21v3 = 21)[name], panel = paste0("cluster ", cluster, " (", table(p$cluster)[as.character(cluster)], " genes)"))
  ggplot(d, aes(months, lfc)) + geom_hline(yintercept = 0, color = "grey60") + geom_line(aes(group = gene_id), color = "grey75", linewidth = 0.3) +
    stat_summary(fun = median, geom = "line", color = "#2a78d6", linewidth = 1.2) + facet_wrap(~ panel, nrow = 1) +
    scale_x_continuous(breaks = c(3, 18, 21)) + labs(x = "age (months)", y = "log2FC vs 3 months")
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

#' DESeq2 datasets behind the four analyses: aging cohort (`exclude` samples dropped), aDR vs AL pair, YMP vs PBS pair.
run_dds <- function(x, samples, exclude = "SRR21355389") {
  grp <- function(g) samples %>% filter(group %in% g) %>% mutate(group = droplevels(group))
  s_age <- samples %>% filter(cohort == "aging", !rownames(.) %in% exclude) %>% mutate(group = droplevels(group))
  list(aging = make_dds(x, s_age, ~ sex + group, min_n = 3),
       aDR = make_dds(x, grp(c("AL", "aDR")), ~ group, min_n = 4),
       YMP = make_dds(x, grp(c("PBS", "YMP")), ~ group, min_n = 3))
}

#' The four result tables from the datasets of run_dds().
de_tables <- function(d, genes) list(
  aging_18v3 = de_table(d$aging, c("group", "aging 18mo", "aging 3mo"), genes),
  aging_21v3 = de_table(d$aging, c("group", "aging 21mo", "aging 3mo"), genes),
  aDR = de_table(d$aDR, c("group", "aDR", "AL"), genes),
  YMP = de_table(d$YMP, c("group", "YMP", "PBS"), genes))

#' STAR GeneCounts matrix saved by notebook 01 (genes x runs, GENCODE M39 versioned gene ids = kallisto's).
load_star_counts <- function(path = file.path(ROOT, "results/star/counts_unstranded.csv.gz")) {
  m <- read_csv(path, show_col_types = FALSE) %>% column_to_rownames(colnames(.)[1]) %>% as.matrix()
  storage.mode(m) <- "integer"; m
}

#' Agreement of two result tables for one contrast (a, b = DE tables of two quantification methods).
compare_contrast <- function(a, b, alpha = 0.05) {
  j <- inner_join(a, b, by = "gene_id", suffix = c("_a", "_b")) %>% filter(!is.na(padj_a), !is.na(padj_b))
  sa <- j$padj_a < alpha; sb <- j$padj_b < alpha; sig <- sa | sb
  tibble(total_DEG_a = sum(a$padj < alpha, na.rm = TRUE), total_DEG_b = sum(b$padj < alpha, na.rm = TRUE),   # all DEGs of each method
         genes_tested_both = nrow(j), DEG_a = sum(sa), DEG_b = sum(sb), shared = sum(sa & sb),   # restricted to genes tested by both
         jaccard = round(sum(sa & sb) / max(sum(sig), 1), 2),
         same_direction = round(mean(sign(j$log2FoldChange_a[sa & sb]) == sign(j$log2FoldChange_b[sa & sb])), 3),
         lfc_spearman_all = round(cor(j$log2FoldChange_a, j$log2FoldChange_b, method = "spearman"), 3),
         lfc_spearman_sig = round(cor(j$log2FoldChange_a[sig], j$log2FoldChange_b[sig], method = "spearman"), 3))
}

#' log2FC of method a (x) against method b (y) per gene; colour = which method calls the gene a DEG.
lfc_scatter <- function(a, b, title, names = c("kallisto", "STAR"), alpha = 0.05) {
  j <- inner_join(a, b, by = "gene_id", suffix = c("_a", "_b")) %>% filter(!is.na(padj_a), !is.na(padj_b)) %>%
    mutate(call = case_when(padj_a < alpha & padj_b < alpha ~ "both", padj_a < alpha ~ paste(names[1], "only"),
                            padj_b < alpha ~ paste(names[2], "only"), TRUE ~ "neither"),
           lab = if_else(call != "neither" & (abs(log2FoldChange_a - log2FoldChange_b) > 1), symbol_a, NA_character_))
  cols <- setNames(c("firebrick", "steelblue", "darkorange", "grey70"), c("both", paste(names[1], "only"), paste(names[2], "only"), "neither"))
  ggplot(j, aes(log2FoldChange_a, log2FoldChange_b, color = call)) + geom_abline(color = "grey50") +
    geom_point(data = filter(j, call == "neither"), size = 0.4, alpha = 0.3) + geom_point(data = filter(j, call != "neither"), size = 1) +
    geom_text_repel(data = filter(j, !is.na(lab)), aes(label = lab), size = 3, show.legend = FALSE, max.overlaps = 20) +
    scale_color_manual(values = cols) + labs(title = title, x = paste("log2FC", names[1]), y = paste("log2FC", names[2]))
}

#' Per sample: Spearman correlation of gene counts for each pair of methods, and total counts (millions) per method.
#' m = named list of gene x sample count matrices on the same genes.
counts_per_sample <- function(m, samples) {
  s <- rownames(samples); pairs <- combn(names(m), 2, simplify = FALSE)
  rho <- set_names(map(pairs, \(p) map_dbl(s, \(i) cor(m[[p[1]]][, i], m[[p[2]]][, i], method = "spearman"))),
                   map_chr(pairs, \(p) paste("spearman", p[1], p[2], sep = "_")))
  tibble(srr = s, group = samples$group, !!!rho, !!!set_names(map(m, \(x) colSums(x[, s]) / 1e6), paste0(names(m), "_M")))
}

#' compare_contrast for every pair of methods and every comparison; de = list(method = list(comparison = DE table)).
method_agreement <- function(de, pairs) {
  bind_rows(lapply(pairs, \(p) bind_rows(lapply(names(de[[p[1]]]), \(n) compare_contrast(de[[p[1]]][[n]], de[[p[2]]][[n]]) %>% mutate(contrast = n))) %>%
    mutate(method_a = p[1], method_b = p[2]))) %>% relocate(method_a, method_b, contrast)
}

#' Genes significant in both / only the first / only the second method of each pair, with their median baseMean.
single_method_calls <- function(de, pairs, alpha = 0.05) {
  one <- \(a, b) inner_join(a, b, by = "gene_id", suffix = c("_a", "_b")) %>% filter(!is.na(padj_a), !is.na(padj_b)) %>%
    mutate(call = case_when(padj_a < alpha & padj_b < alpha ~ "both", padj_a < alpha ~ "first only", padj_b < alpha ~ "second only", TRUE ~ "neither")) %>%
    filter(call != "neither") %>% group_by(call) %>% summarise(n_genes = n(), median_baseMean = round(median((baseMean_a + baseMean_b) / 2)), .groups = "drop")
  bind_rows(lapply(pairs, \(p) bind_rows(lapply(names(de[[p[1]]]), \(n) one(de[[p[1]]][[n]], de[[p[2]]][[n]]) %>% mutate(contrast = n))) %>%
    mutate(first = p[1], second = p[2]))) %>% relocate(first, second, contrast)
}

#' One row per method with the criteria of docs/pipeline.md ("Choosing the counting method"). counts = list of tximport lists / count
#' matrices; dds, de = run_dds / de_tables per method; ids = genes present in all methods.
method_evidence <- function(counts, dds, de, samples, genes, ids) {
  s <- rownames(samples)
  # usable data: share of the trimmed read pairs counted in a gene; genes with >= 10 counts in >= 3 runs
  inp <- read_csv(file.path(ROOT, "results/star/star_qc.csv"), show_col_types = FALSE) %>% column_to_rownames("SRR")
  pairs_in <- inp[s, "input_pairs_M"] * 1e6
  cnt <- lapply(counts, \(m) if (is.matrix(m)) m[, s] else m$counts[, s])
  usable <- bind_rows(lapply(names(cnt), \(m) tibble(method = m, pct = 100 * colSums(cnt[[m]]) / pairs_in))) %>% group_by(method) %>%
    summarise(counted_pairs_pct_median = round(median(pct), 1), min = round(min(pct), 1), max = round(max(pct), 1)) %>%
    left_join(tibble(method = names(cnt), genes_detected = map_int(cnt, \(m) sum(rowSums(m >= 10) >= 3))), by = "method")
  # replicate agreement: median over groups of the median pairwise Pearson r of log2 CPM (genes with mean >= 10 counts in every method)
  keep <- ids[Reduce(`&`, lapply(cnt, \(m) rowMeans(m[ids, ]) >= 10))]
  lcpm <- lapply(cnt, \(m) { m <- m[keep, ]; log2(sweep(m, 2, colSums(m) / 1e6, "/") + 1) })
  repl <- bind_rows(lapply(names(lcpm), \(m) tibble(method = m, r = map_dbl(split(s, samples$group), \(g) { cm <- cor(lcpm[[m]][, g]); median(cm[lower.tri(cm)]) })))) %>%
    group_by(method) %>% summarise(replicate_r_median = round(median(r), 4))
  # noise: median DESeq2 dispersion per analysis over genes tested by all methods
  tested <- Reduce(intersect, lapply(dds, \(d) Reduce(intersect, lapply(d, rownames))))
  disp <- bind_rows(lapply(names(dds), \(m) tibble(method = m, !!!lapply(dds[[m]], \(d) round(median(setNames(dispersions(d), rownames(d))[intersect(rownames(d), tested)], na.rm = TRUE), 3)))))
  # stability: aging DEGs (21 vs 3 months) without / with the outlier SRR21355389
  s_all <- samples %>% filter(cohort == "aging") %>% mutate(group = droplevels(group))
  stab <- bind_rows(lapply(names(counts), \(m) {
    all12 <- de_table(make_dds(counts[[m]], s_all, ~ sex + group, min_n = 4), c("group", "aging 21mo", "aging 3mo"), genes)
    cd <- compare_de(de[[m]]$aging_21v3, all12)
    tibble(method = m, DEG_without = cd$DEG_a, DEG_with = cd$DEG_b, jaccard_outlier = round(cd$DEG_both / (cd$DEG_a + cd$DEG_b - cd$DEG_both), 2)) }))
  # positive controls (the paper: 5 inflammation genes up with aging, 23 circadian/cycling genes up with aDR)
  up_aging <- c("C4b", "Trem2", "Ctss", "Tyrobp", "Gfap")
  adr23 <- c("Chordc1", "Cxcl12", "Fabp7", "Kdr", "Adipor2", "Gjb6", "Mertk", "Cry1", "Klf15", "Dbp", "Usp2", "Etnppl", "Ciart", "Sult1a1",
             "Zbtb16", "Fkbp5", "Sgk1", "Tsc22d3", "Cirbp", "Rbm3", "Plin4", "Fmo2", "Hif3a")
  ctrl <- bind_rows(lapply(names(de), \(m) tibble(method = m,
    aging_genes_up_of_5 = sum(filter(de[[m]]$aging_21v3, symbol %in% up_aging)$log2FoldChange > 0),
    aDR_genes_up_of_23 = sum(filter(de[[m]]$aDR, symbol %in% adr23)$log2FoldChange > 0))))
  usable %>% left_join(repl, by = "method") %>% left_join(rename_with(disp, ~ paste0("dispersion_", .x), -method), by = "method") %>%
    left_join(stab, by = "method") %>% left_join(ctrl, by = "method")
}

#' Technical vs biological variation on TPM (TA's design). a, b = gene x sample TPM of two methods (same genes and samples).
#' Genes: mean TPM >= min_tpm in both methods. All values are log2((TPM + 1) / (TPM + 1)) per gene, centred on their median
#' (`shift`): TPMs sum to a million, so a few abundant genes shift all others (removed by DESeq2's normalisation as well):
#'   technical  = method a vs b on the same mouse; "offset removed" subtracts each gene's mean method difference over all mice
#'                (a fixed per-gene offset cancels when both mice are counted with the same method);
#'   biological = two mice of the same group and sex, one method; pair_lfc = the two methods' log2FC for the same pair of mice, subtracted.
#' Returns per-comparison spreads (`per`) and the pooled log2FC values for plotting (`lfc`).
tpm_variation <- function(a, b, samples, names = c("kallisto", "RSEM"), min_tpm = 1, exclude = "SRR21355389") {
  keep <- rowMeans(a) >= min_tpm & rowMeans(b) >= min_tpm
  la <- log2(a[keep, ] + 1); lb <- log2(b[keep, ] + 1)
  d <- la - lb; d0 <- d - rowMeans(d)
  level <- cut(log10((rowMeans(a[keep, ]) + rowMeans(b[keep, ])) / 2), c(-Inf, 1, 2, Inf), labels = c("TPM 1-10", "TPM 10-100", "TPM > 100"))
  pairs <- samples %>% rownames_to_column("srr") %>% filter(!srr %in% exclude) %>% group_by(group, sex) %>% filter(n() >= 2) %>%
    summarise(p = list(combn(srr, 2, simplify = FALSE)), .groups = "drop") %>% pull(p) %>% unlist(recursive = FALSE)
  kinds <- c(setNames(list(lapply(colnames(d), \(s) d[, s]), lapply(colnames(d), \(s) d0[, s])),
                      paste(c("technical:", "technical, offset removed:"), names[1], "vs", names[2], "(same mouse)")),
             setNames(lapply(list(lb, la), \(l) lapply(pairs, \(p) l[, p[1]] - l[, p[2]])),
                      paste0("biological: two mice, same group and sex (", names[2:1], ")")),
             setNames(list(lapply(pairs, \(p) (la[, p[1]] - la[, p[2]]) - (lb[, p[1]] - lb[, p[2]]))),
                      paste("pair log2FC,", names[1], "minus", names[2])))
  kinds <- map(kinds, \(xs) map(xs, \(x) structure(x - median(x), shift = median(x))))
  per <- imap_dfr(kinds, \(xs, k) map_dfr(xs, \(x) tibble(comparison = k, shift = attr(x, "shift"), sd = sd(x), median_abs = median(abs(x)),
                                                         over_2fold = mean(abs(x) > 1), !!!map(split(abs(x), level), median))))
  list(per = per, genes = sum(keep),
       lfc = imap_dfr(kinds[c(1, 2, 3)], \(xs, k) tibble(comparison = k, log2FC = unlist(xs))))
}

#' Summary of tpm_variation: median over comparisons of each spread (log2 units; over_2fold in % of genes).
tpm_variation_table <- function(v) v$per %>% group_by(comparison = fct_inorder(comparison)) %>%
  summarise(n = n(), max_abs_shift = round(max(abs(shift)), 2), across(c(sd, median_abs, starts_with("TPM")), \(x) round(median(x), 3)), over_2fold_pct = round(100 * median(over_2fold), 1))

#' Pooled distributions of log2FC: technical (with and without the per-gene offset) vs biological; ideally all centred on 0
#' and the technical narrower.
tpm_variation_plot <- function(v, xlim = c(-2, 2)) {
  cols <- setNames(c("#2a78d6", "#2a78d6", "#eb6834"), unique(v$lfc$comparison)); lty <- setNames(c("dotted", "solid", "solid"), names(cols))
  ggplot(v$lfc, aes(log2FC, color = comparison, linetype = comparison)) + geom_vline(xintercept = 0, color = "grey60") +
    geom_density(linewidth = 1, adjust = 1.5) + coord_cartesian(xlim = xlim) +
    scale_color_manual(values = cols) + scale_linetype_manual(values = lty) +
    labs(x = "log2 fold change per gene", y = "density", color = NULL, linetype = NULL) +
    theme(legend.position = "bottom", legend.direction = "vertical")
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
