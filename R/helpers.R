# R helpers for notebooks 03, 04, 05 (run them from notebooks/; paths are relative to the repo root).
#
# Contents (functions grouped by purpose)
#   Load data ........ load_samples, load_genes, load_rsem, load_txi (kallisto), subset_txi
#   DESeq2 ........... make_dds, de_table, n_deg, run_dds, de_tables
#   Explore (03) ..... COL, pca_plot, outlier_genes, outlier_plot, volcano_plot, age_lfc_plot, deg_count_plot, annotate_degs
#   Compare methods .. compare_de, compare_contrast, lfc_scatter, method_overlap_plot   (any two methods; outlier sensitivity)
#                      counts_per_sample, method_agreement, single_method_calls, method_evidence, polyg_effect (notebook 04)
#                      tpm_variation, tpm_variation_table, tpm_variation_plot (technical vs biological variation)
#   Aging vs treatment overlap_with_aging, signature_score, score_test, score_plot, lfc_vs_aging_plot,
#                      msigdb_sets, run_gsea, nes_tiles, gsea_concordance, gene_table, shared_with_aging
#   Team common setup  lfc_of, reversal_measure, write_common (results/common/)
# Terms: docs/glossary.md. Conventions: `res` = DESeq2 results table (columns gene_id, baseMean, log2FoldChange, lfcSE, stat,
# pvalue, padj, symbol); `x` = tximport list (kallisto, RSEM) or gene x sample count matrix (STAR); `dds` = DESeq2 dataset.
suppressPackageStartupMessages({library(tidyverse); library(tximport); library(DESeq2); library(ggrepel)})
ROOT <- normalizePath("..")
# One look for every figure: very large text, no minor grid, tight margins (notebooks set the figure size per plot).
theme_set(theme_minimal(base_size = 40) +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "#333333"),
        plot.title = element_text(size = rel(1)), plot.subtitle = element_text(size = rel(0.85), colour = "#444444"),
        strip.text = element_text(size = rel(0.95)), legend.title = element_text(size = rel(0.9)),
        legend.key.size = unit(2.2, "lines"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0.4, "lines"),
        panel.spacing = unit(2, "lines"), plot.margin = margin(4, 8, 4, 4)))
LAB <- 13   # size of text labels inside plots (mm; ~37 pt)

#' The 27 runs with cohort / experiment (aging, aDR, YMP) / sex / age / group (factors ordered with the control first); age_class = young (3 months) / old (18 + 21 months)
#' in the aging cohort, NA otherwise.
load_samples <- function() {
  read_csv(file.path(ROOT, "metadata/samples.csv"), show_col_types = FALSE) %>%
    transmute(srr = SRR, mouse_id, sex, age = age_months, group,
              cohort = if_else(BioProject == "PRJNA875066", "aging", "rejuvenation")) %>%
    mutate(group = factor(group, levels = c("aging 3mo", "aging 18mo", "aging 21mo", "AL", "aDR", "PBS", "YMP")),
           sex = factor(sex),
           age_class = factor(case_when(group == "aging 3mo" ~ "young", group %in% c("aging 18mo", "aging 21mo") ~ "old"), levels = c("young", "old")),
           experiment = factor(case_when(cohort == "aging" ~ "aging", group %in% c("AL", "aDR") ~ "aDR", TRUE ~ "YMP"), levels = c("aging", "aDR", "YMP"))) %>%
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

#' Tidy results for one contrast (c(factor, treated, control); control = the reference level), annotated with gene symbols.
#' alpha = 0.05 = our DEG cutoff, so DESeq2's independent filtering is tuned for padj < 0.05 (team setup).
#' shrink = TRUE adds lfc_shrunk (apeglm), used for plots and for comparing fold changes; p-values stay unshrunken.
de_table <- function(dds, contrast, genes, shrink = FALSE) {
  r <- results(dds, contrast = contrast, alpha = 0.05) %>% as.data.frame() %>% rownames_to_column("gene_id")
  if (shrink) {
    coef <- make.names(paste(contrast[1], contrast[2], "vs", contrast[3], sep = "_"))
    s <- lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE)
    r$lfc_shrunk <- s[r$gene_id, "log2FoldChange"]
  }
  r %>% left_join(genes, by = "gene_id") %>% arrange(padj)
}

#' Fold change to plot or compare: apeglm-shrunken when the table has it (team setup), else DESeq2's estimate.
lfc_of <- function(r) if ("lfc_shrunk" %in% names(r)) r$lfc_shrunk else r$log2FoldChange

n_deg <- function(res, alpha = 0.05) res %>% filter(padj < alpha) %>% summarise(up = sum(log2FoldChange > 0), down = sum(log2FoldChange < 0))

# Colours (validated palette): experiment = categorical slots 1-3; age = one blue ramp (light = young); up / down = red / blue poles.
COL <- list(experiment = c(aging = "#2a78d6", aDR = "#eb6834", YMP = "#1baf7a"),
            age = c(`aging 3mo` = "#86b6ef", `aging 18mo` = "#2a78d6", `aging 21mo` = "#104281"),
            dir = c(up = "#e34948", down = "#2a78d6", `not DEG` = "#cfcfcf"))

#' PCA of the variance-stabilised counts (vst, blind to the groups). `label` = runs to name on the plot (default none).
pca_plot <- function(dds, color, shape = NULL, label = character(0), colors = NULL) {
  v <- suppressMessages(vst(dds, blind = TRUE))
  d <- suppressMessages(plotPCA(v, intgroup = color, returnData = TRUE)); pv <- round(100 * attr(d, "percentVar"))
  d[[color]] <- colData(v)[[color]]
  if (!is.null(shape)) d[[shape]] <- colData(v)[[shape]]
  ggplot(d, aes(PC1, PC2, color = .data[[color]])) +
    (if (is.null(shape)) geom_point(size = 8) else geom_point(aes(shape = .data[[shape]]), size = 8)) +
    geom_text_repel(data = filter(d, name %in% label), aes(label = name), size = LAB, color = "black",
                    box.padding = 1.2, seed = 42, show.legend = FALSE) +
    (if (!is.null(colors)) scale_color_manual(values = colors) else NULL) +
    labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"), color = NULL, shape = NULL)
}

#' Genes that most distinguish one sample from the rest (vst difference ≈ log2 difference) - to diagnose an outlier.
outlier_genes <- function(dds, sample, genes, n = 10) {
  v <- assay(suppressMessages(vst(dds, blind = TRUE)))
  tibble(gene_id = rownames(v), vst_diff = v[, sample] - rowMeans(v[, colnames(v) != sample])) %>%
    left_join(genes, by = "gene_id") %>% arrange(desc(abs(vst_diff))) %>% head(n) %>% select(symbol, vst_diff)
}

#' Bars: how many times higher each of the outlier's top genes is than in the other mice (2^vst_diff).
outlier_plot <- function(og, sample) {
  ggplot(og, aes(2^vst_diff, fct_reorder(symbol, vst_diff))) + geom_col(fill = "#2a78d6", width = 0.7) +
    geom_text(aes(label = paste0(round(2^vst_diff), "×")), hjust = -0.15, size = LAB) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = paste("× higher in", sample, "than in other aging mice"), y = NULL) +
    theme(panel.grid.major.y = element_blank())
}

#' Volcano plot(s): log2FC vs -log10 padj; red / blue = DEGs up / down (padj < alpha); the n_label most significant DEGs are named.
#' tables = named list of DE tables (one panel each); the panel title gives the DEG counts.
volcano_plot <- function(tables, alpha = 0.05, n_label = 8) {
  d <- bind_rows(lapply(names(tables), \(n) {
    r <- filter(tables[[n]], !is.na(padj)); k <- n_deg(r, alpha)
    mutate(r, lfc = lfc_of(r), panel = sprintf("%s: %d up, %d down", n, k$up, k$down),
           call = factor(case_when(padj >= alpha ~ "not DEG", log2FoldChange > 0 ~ "up", TRUE ~ "down"), names(COL$dir)))
  })) %>% mutate(panel = fct_inorder(panel))
  lab <- d %>% filter(call != "not DEG") %>% group_by(panel) %>% slice_min(padj, n = n_label, with_ties = FALSE)
  ggplot(d, aes(lfc, -log10(padj), color = call)) +
    geom_point(data = filter(d, call == "not DEG"), size = 1.2) + geom_point(data = filter(d, call != "not DEG"), size = 3.5) +
    geom_text_repel(data = lab, aes(label = symbol), size = LAB, color = "black", bg.color = "white", bg.r = 0.12, seed = 42, max.overlaps = Inf, min.segment.length = 0) +
    scale_color_manual(values = COL$dir, drop = FALSE) + facet_wrap(~ panel, scales = "free") +
    labs(x = "log2 fold change (shrunken)", y = "-log10 adjusted p", color = NULL) +
    theme(legend.position = "bottom") + guides(color = guide_legend(override.aes = list(size = 7)))
}

#' Aging DEGs: log2FC at 18 months vs at 21 months (each vs 3 months); points on the diagonal changed equally by both ages.
#' age3 = list(`18v3`, `21v3`) of DE tables; degs = the aging DEGs (rows of the old vs young table).
age_lfc_plot <- function(age3, degs, n_label = 6) {
  lfc <- \(t) setNames(lfc_of(t), t$gene_id)[degs$gene_id]
  d <- degs %>% transmute(gene_id, symbol, call = factor(if_else(log2FoldChange > 0, "up", "down"), c("up", "down")),
                          lfc18 = lfc(age3$`18v3`), lfc21 = lfc(age3$`21v3`))
  lab <- slice_max(d, abs(lfc21), n = n_label)
  ggplot(d, aes(lfc18, lfc21, color = call)) + geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
    geom_abline(linetype = "dashed", color = "grey50") + geom_point(size = 4) +
    geom_text_repel(data = lab, aes(label = symbol), size = LAB, color = "black", bg.color = "white", bg.r = 0.12, seed = 42) +
    scale_color_manual(values = COL$dir) +
    labs(x = "log2FC 18 vs 3 months", y = "log2FC 21 vs 3 months", color = "old vs young",
         subtitle = sprintf("%d aging DEGs · Pearson r = %.2f", nrow(d), cor(d$lfc18, d$lfc21))) +
    theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1))
}

#' DEG counts per comparison: up to the right, down to the left. tables = named list of DE tables (names used as labels).
deg_count_plot <- function(tables, alpha = 0.05) {
  d <- bind_rows(lapply(tables, n_deg, alpha = alpha), .id = "comparison") %>% mutate(comparison = fct_rev(fct_inorder(comparison))) %>%
    pivot_longer(c(up, down), names_to = "call", values_to = "n") %>% mutate(x = if_else(call == "up", n, -n))
  ggplot(d, aes(x, comparison, fill = call)) + geom_col(width = 0.6) + geom_vline(xintercept = 0, color = "grey40") +
    geom_text(aes(label = n, hjust = if_else(call == "up", -0.3, 1.3)), size = LAB) +
    scale_fill_manual(values = COL$dir[c("up", "down")]) +
    scale_x_continuous(labels = abs, expand = expansion(mult = 0.2)) +
    labs(x = "DEGs (padj < 0.05):  ← down · up →", y = NULL, fill = NULL) +
    theme(panel.grid.major.y = element_blank(), legend.position = "none")
}

#' Sensitivity of one contrast to a choice: DEG counts, overlap and LFC agreement between two result tables.
compare_de <- function(a, b, alpha = 0.05) {
  j <- inner_join(a, b, by = "gene_id", suffix = c("_a", "_b")) %>% filter(!is.na(padj_a), !is.na(padj_b))
  sig <- j$padj_a < alpha | j$padj_b < alpha
  tibble(DEG_a = sum(j$padj_a < alpha), DEG_b = sum(j$padj_b < alpha), DEG_both = sum(j$padj_a < alpha & j$padj_b < alpha),
         lfc_spearman = cor(j$log2FoldChange_a[sig], j$log2FoldChange_b[sig], method = "spearman"))
}

#' DESeq2 datasets behind the three analyses: aging cohort old vs young (`exclude` samples dropped), aDR vs AL pair, YMP vs PBS pair.
run_dds <- function(x, samples, exclude = "SRR21355389") {
  grp <- function(g) samples %>% filter(group %in% g) %>% mutate(group = droplevels(group))
  s_age <- samples %>% filter(cohort == "aging", !rownames(.) %in% exclude) %>% mutate(group = droplevels(group))
  list(aging = make_dds(x, s_age, ~ sex + age_class, min_n = 4),
       aDR = make_dds(x, grp(c("AL", "aDR")), ~ group, min_n = 4),
       YMP = make_dds(x, grp(c("PBS", "YMP")), ~ group, min_n = 3))
}

#' The three result tables from the datasets of run_dds(); aging = old (18 + 21 months) vs young (3 months).
de_tables <- function(d, genes, shrink = FALSE) list(
  aging = de_table(d$aging, c("age_class", "old", "young"), genes, shrink),
  aDR = de_table(d$aDR, c("group", "aDR", "AL"), genes, shrink),
  YMP = de_table(d$YMP, c("group", "YMP", "PBS"), genes, shrink))

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

#' log2FC of method a (x) against method b (y), one panel per comparison; colour = which method calls the gene a DEG.
#' a, b = named lists of DE tables (same names).
lfc_scatter <- function(a, b, names = c("kallisto", "RSEM"), alpha = 0.05) {
  lv <- c("both", paste(names, "only"), "neither")
  d <- bind_rows(lapply(names(a), \(n) inner_join(a[[n]], b[[n]], by = "gene_id", suffix = c("_a", "_b")) %>%
    filter(!is.na(padj_a), !is.na(padj_b)) %>% mutate(panel = n))) %>% mutate(panel = fct_inorder(panel),
    call = factor(case_when(padj_a < alpha & padj_b < alpha ~ lv[1], padj_a < alpha ~ lv[2], padj_b < alpha ~ lv[3], TRUE ~ lv[4]), lv))
  ggplot(d, aes(log2FoldChange_a, log2FoldChange_b, color = call)) + geom_abline(color = "grey50") +
    geom_point(data = filter(d, call == "neither"), size = 0.8) + geom_point(data = filter(d, call != "neither"), size = 4) +
    scale_color_manual(values = setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#d6d6d6"), lv), breaks = lv[1:3]) + facet_wrap(~ panel) +
    coord_cartesian(xlim = c(-4, 4), ylim = c(-4, 4)) +
    labs(x = paste("log2FC", names[1]), y = paste("log2FC", names[2]), color = NULL) +
    theme(legend.position = "bottom") + guides(color = guide_legend(override.aes = list(size = 7)))
}

#' DEGs per comparison split into called by both methods / only one (stacked bars, counts written in the bars).
#' cmp = method_agreement() rows (columns contrast, DEG_a, DEG_b, shared).
method_overlap_plot <- function(cmp, names = c("kallisto", "RSEM")) {
  lv <- c("both", paste(names, "only"))
  d <- cmp %>% transmute(contrast = fct_rev(fct_inorder(contrast)), both = shared, a = DEG_a - shared, b = DEG_b - shared) %>%
    pivot_longer(-contrast, names_to = "call", values_to = "n") %>% mutate(call = factor(setNames(lv, c("both", "a", "b"))[call], lv))
  ggplot(d, aes(n, contrast, fill = call)) + geom_col(width = 0.7, color = "white", linewidth = 1, position = position_stack(reverse = TRUE)) +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5, reverse = TRUE), size = LAB, color = "white", fontface = "bold") +
    scale_fill_manual(values = setNames(c("#2a78d6", "#eb6834", "#1baf7a"), lv)) + scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = "DEGs (padj < 0.05, genes tested by both)", y = NULL, fill = NULL) +
    theme(legend.position = "top", panel.grid.major.y = element_blank())
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
  # stability: aging DEGs (old vs young) without / with the outlier SRR21355389
  s_all <- samples %>% filter(cohort == "aging") %>% mutate(group = droplevels(group))
  stab <- bind_rows(lapply(names(counts), \(m) {
    all12 <- de_table(make_dds(counts[[m]], s_all, ~ sex + age_class, min_n = 4), c("age_class", "old", "young"), genes)
    cd <- compare_de(de[[m]]$aging, all12)
    tibble(method = m, DEG_without = cd$DEG_a, DEG_with = cd$DEG_b, jaccard_outlier = round(cd$DEG_both / (cd$DEG_a + cd$DEG_b - cd$DEG_both), 2)) }))
  # positive controls (the paper: 5 inflammation genes up with aging, 23 circadian/cycling genes up with aDR)
  up_aging <- c("C4b", "Trem2", "Ctss", "Tyrobp", "Gfap")
  adr23 <- c("Chordc1", "Cxcl12", "Fabp7", "Kdr", "Adipor2", "Gjb6", "Mertk", "Cry1", "Klf15", "Dbp", "Usp2", "Etnppl", "Ciart", "Sult1a1",
             "Zbtb16", "Fkbp5", "Sgk1", "Tsc22d3", "Cirbp", "Rbm3", "Plin4", "Fmo2", "Hif3a")
  ctrl <- bind_rows(lapply(names(de), \(m) tibble(method = m,
    aging_genes_up_of_5 = sum(filter(de[[m]]$aging, symbol %in% up_aging)$log2FoldChange > 0),
    aDR_genes_up_of_23 = sum(filter(de[[m]]$aDR, symbol %in% adr23)$log2FoldChange > 0))))
  usable %>% left_join(repl, by = "method") %>% left_join(rename_with(disp, ~ paste0("dispersion_", .x), -method), by = "method") %>%
    left_join(stab, by = "method") %>% left_join(ctrl, by = "method")
}

#' Effect of poly-G trimming on counts. with / without = gene x sample counts with / without poly-G trimming (same samples).
#' Genes with mean count >= min_mean are compared. Returns per sample (`per`): Spearman r, total counts (millions), % of genes changing
#' > 5 % and > 2-fold; and per gene (`top`): median log2(without / with) over the samples, largest changes first.
polyg_effect <- function(with, without, genes, min_mean = 10) {
  ids <- intersect(rownames(with), rownames(without)); a <- with[ids, , drop = FALSE]; b <- without[ids, colnames(a), drop = FALSE]
  keep <- (rowMeans(a) + rowMeans(b)) / 2 >= min_mean
  l <- log2((b[keep, , drop = FALSE] + 1) / (a[keep, , drop = FALSE] + 1))
  list(genes_compared = sum(keep),
       per = tibble(srr = colnames(a), spearman = round(map_dbl(colnames(a), \(i) cor(a[keep, i], b[keep, i], method = "spearman")), 5),
                    counts_with_M = round(colSums(a) / 1e6, 3), counts_without_M = round(colSums(b) / 1e6, 3),
                    genes_changed_5pct = round(100 * colMeans(abs(l) > log2(1.05)), 3), genes_changed_2fold = round(100 * colMeans(abs(l) > 1), 3)),
       top = tibble(gene_id = rownames(l), log2_without_vs_with = round(apply(l, 1, median), 2),
                    mean_with = round(rowMeans(a[keep, , drop = FALSE])), mean_without = round(rowMeans(b[keep, , drop = FALSE]))) %>%
         left_join(genes, by = "gene_id") %>% relocate(symbol, .after = gene_id) %>% arrange(desc(abs(log2_without_vs_with))))
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

#' Pooled distributions of log2FC: two mice (biological) vs the same mouse counted by both methods (technical, raw and offset removed).
tpm_variation_plot <- function(v, xlim = c(-1, 1)) {
  k <- unique(v$lfc$comparison)   # technical, technical offset removed, biological (RSEM)
  lab <- setNames(c("same mouse, two methods", "same mouse, two methods, offset removed", "two mice, one method (RSEM)"), k)
  d <- mutate(v$lfc, comparison = factor(lab[comparison], lab[c(3, 2, 1)]))
  ggplot(d, aes(log2FC, color = comparison, linetype = comparison)) + geom_vline(xintercept = 0, color = "grey60") +
    geom_density(linewidth = 2.5, adjust = 1.5) + coord_cartesian(xlim = xlim) +
    scale_color_manual(values = setNames(c("#eb6834", "#2a78d6", "#2a78d6"), lab[c(3, 2, 1)])) +
    scale_linetype_manual(values = setNames(c("solid", "solid", "dotted"), lab[c(3, 2, 1)])) +
    labs(x = "log2 fold change per gene", y = "density", color = NULL, linetype = NULL) +
    theme(legend.position = "bottom", legend.direction = "vertical", legend.key.width = unit(4, "lines"))
}

# ---- aging vs rejuvenation ------------------------------------------------------------------

#' All DEGs (padj < alpha) of a named list of DE tables, one row per gene and comparison, with the gene's full name and type from
#' org.Mm.eg.db (NCBI Gene; NA = no entry, mostly predicted genes and lncRNAs).
annotate_degs <- function(tables, alpha = 0.05) {
  d <- bind_rows(lapply(names(tables), \(n) filter(tables[[n]], padj < alpha) %>% mutate(comparison = n)))
  db <- suppressPackageStartupMessages(org.Mm.eg.db::org.Mm.eg.db)
  q <- \(keys, type) suppressMessages(AnnotationDbi::select(db, keys, c("SYMBOL", "GENENAME", "GENETYPE"), type))
  # by Ensembl id, kept only where NCBI's symbol matches ours (ids can map to another gene between releases); else by symbol
  by_id <- q(unique(sub("\\..*", "", d$gene_id)), "ENSEMBL") %>% filter(!is.na(SYMBOL))
  by_sym <- q(unique(d$symbol), "SYMBOL") %>% filter(!is.na(GENENAME))
  ann <- distinct(d, gene_id, symbol) %>% mutate(ENSEMBL = sub("\\..*", "", gene_id)) %>%
    left_join(semi_join(by_id, ., by = c("ENSEMBL", "SYMBOL" = "symbol")) %>% distinct(ENSEMBL, .keep_all = TRUE) %>% select(ENSEMBL, n1 = GENENAME, t1 = GENETYPE), by = "ENSEMBL") %>%
    left_join(distinct(by_sym, SYMBOL, .keep_all = TRUE) %>% select(symbol = SYMBOL, n2 = GENENAME, t2 = GENETYPE), by = "symbol") %>%
    transmute(gene_id, name = coalesce(n1, n2), type = coalesce(t1, t2))
  d %>% left_join(ann, by = "gene_id") %>%
    transmute(comparison, symbol, name, type, log2FC = round(log2FoldChange, 2), padj = signif(padj, 2), direction = if_else(log2FoldChange > 0, "up", "down"), gene_id) %>%
    arrange(factor(comparison, names(tables)), padj)
}

#' Treatment DEGs that are also aging DEGs (old vs young): log2FC in both; reversal = opposite sign to aging.
shared_with_aging <- function(res, treatments = c("aDR", "YMP"), alpha = 0.05) {
  aging <- filter(res$aging, padj < alpha) %>% select(gene_id, symbol, lfc_aging = log2FoldChange)
  bind_rows(lapply(treatments, \(n) inner_join(aging, filter(res[[n]], padj < alpha) %>% select(gene_id, lfc_treatment = log2FoldChange), by = "gene_id") %>%
    mutate(treatment = n))) %>% mutate(across(starts_with("lfc"), \(x) round(x, 2)), vs_aging = if_else(sign(lfc_aging) == sign(lfc_treatment), "same", "reversal")) %>%
    select(treatment, symbol, lfc_aging, lfc_treatment, vs_aging)
}
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

#' MSigDB mouse-native gene sets (db_species = "MM") as named lists of gene symbols: Hallmark (MH, 50) and GO biological process (M5).
#' The msigdbr version is printed so it is recorded in the notebook output.
msigdb_sets <- function() {
  to_list <- \(d) split(d$gene_symbol, d$gs_name)
  cat("msigdbr", as.character(packageVersion("msigdbr")), "| source: mouse-native MSigDB (db_species = MM)\n")
  suppressMessages(list(hallmark = to_list(msigdbr::msigdbr(species = "Mus musculus", db_species = "MM", collection = "MH")),
                        gobp = to_list(msigdbr::msigdbr(species = "Mus musculus", db_species = "MM", collection = "M5", subcollection = "GO:BP"))))
}

#' fgsea on the DESeq2 Wald statistic (one value per gene symbol).
run_gsea <- function(res, sets, seed = 42) {   # seed 42 = team setup
  st <- res %>% filter(!is.na(stat), !is.na(symbol), symbol != "") %>% group_by(symbol) %>% slice_max(abs(stat), n = 1, with_ties = FALSE) %>% ungroup()
  set.seed(seed)
  withCallingHandlers(fgsea::fgsea(sets, setNames(st$stat, st$symbol), minSize = 15, maxSize = 500),
                      warning = \(w) if (grepl("ties in the preranked", conditionMessage(w))) invokeRestart("muffleWarning")) %>%   # ~0.1 % tied stats: harmless
    as_tibble() %>% select(-leadingEdge)
}

#' Treatment log2FC against aging log2FC (genes tested in both), one panel per treatment; genes that are DEGs of both are named.
#' Upper-left / lower-right = opposite to aging (reversal). res = list(aging, aDR, YMP) of DE tables.
lfc_vs_aging_plot <- function(res, treatments = c("aDR", "YMP"), alpha = 0.05) {
  lv <- c("both", "aging only", "treatment only", "neither")
  d <- bind_rows(lapply(treatments, \(n) inner_join(res$aging, res[[n]], by = "gene_id", suffix = c("_age", "_int")) %>%
    filter(!is.na(padj_age), !is.na(padj_int)) %>% mutate(panel = n, x = lfc_of(res$aging)[match(gene_id, res$aging$gene_id)],
                                                         y = lfc_of(res[[n]])[match(gene_id, res[[n]]$gene_id)]))) %>%
    mutate(call = factor(case_when(padj_age < alpha & padj_int < alpha ~ "both", padj_age < alpha ~ "aging only",
                                   padj_int < alpha ~ "treatment only", TRUE ~ "neither"), lv))
  ggplot(d, aes(x, y, color = call)) +
    geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
    geom_point(data = filter(d, call == "neither"), size = 0.8) + geom_point(data = filter(d, call %in% lv[2:3]), size = 3.5) +
    geom_point(data = filter(d, call == "both"), size = 7) +
    geom_text_repel(data = filter(d, call == "both"), aes(label = symbol_age), size = LAB, color = "black", bg.color = "white", bg.r = 0.12, seed = 42, box.padding = 1, min.segment.length = 0) +
    scale_color_manual(values = setNames(c("#e34948", "#2a78d6", "#eb6834", "#d6d6d6"), lv), breaks = lv[1:3]) + facet_wrap(~ panel) +
    coord_cartesian(xlim = c(-2, 3), ylim = c(-2, 2)) +
    labs(x = "log2FC aging (old vs young, shrunken)", y = "log2FC treatment (shrunken)", color = NULL) +
    theme(legend.position = "bottom") + guides(color = guide_legend(override.aes = list(size = 7)))
}

#' Per-mouse signature scores, one panel per experiment (aging, aDR, YMP); `tests` = score_test() rows, shown in the panel titles.
score_plot <- function(sc, tests) {
  ttl <- c(aging = "aging (source)", setNames(sprintf("%s: %+.2f, p = %.2g", tests$treated, tests$diff, tests$t_test_p), tests$treated))
  d <- bind_rows(aging = sc$aging, aDR = sc$aDR, YMP = sc$YMP, .id = "panel") %>%
    mutate(group = factor(sub("aging ", "", group), c("3mo", "18mo", "21mo", "AL", "aDR", "PBS", "YMP")),
           panel = factor(ttl[panel], ttl[c("aging", "aDR", "YMP")]), exp = names(ttl)[match(panel, ttl)])
  ggplot(d, aes(group, score, color = exp)) + geom_hline(yintercept = 0, color = "grey70") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 1, color = "grey30") +
    geom_jitter(width = 0.12, height = 0, size = 8) + scale_color_manual(values = COL$experiment) +
    facet_wrap(~ panel, scales = "free_x") + labs(x = NULL, y = "aging score (higher = older)") +
    theme(legend.position = "none", panel.grid.major.x = element_blank())
}

#' Hallmark NES tiles (blue = down, red = up; * = padj < alpha) for the pathways significant in aging, ordered by the aging NES.
nes_tiles <- function(g, alpha = 0.05) {
  d <- bind_rows(g, .id = "contrast") %>% mutate(pathway = str_to_sentence(gsub("_", " ", sub("^HALLMARK_", "", pathway))) %>%
    str_replace_all(c("^Kras" = "KRAS", "^Tnfa signaling via nfkb" = "TNFα signaling via NF-κB", "^Mtorc1" = "mTORC1", "^Myc" = "MYC", "^Il2 stat5" = "IL2-STAT5",
                      "^Il6 jak stat3" = "IL6-JAK-STAT3", "^E2f" = "E2F", "^P53" = "p53", "^Uv" = "UV", "^Dna" = "DNA", "^Tgf beta" = "TGF-β", "^Pi3k akt mtor" = "PI3K-AKT-mTOR",
                      "^Wnt beta catenin" = "WNT-β-catenin", "^G2m" = "G2M", "^Mtorc1" = "mTORC1")))
  keep <- d$pathway[d$contrast == "aging" & d$padj < alpha]; ord <- d %>% filter(contrast == "aging") %>% arrange(NES) %>% pull(pathway)
  d <- d %>% filter(pathway %in% keep) %>% mutate(pathway = factor(pathway, ord), contrast = factor(contrast, c("aging", "aDR", "YMP")))
  ggplot(d, aes(contrast, pathway, fill = NES)) + geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = ifelse(padj < alpha, "*", "")), size = 18, vjust = 0.75) +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948") + scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) + theme(panel.grid = element_blank(), legend.key.height = unit(3, "lines"))
}

#' Pathways significant in aging: how many does each intervention reverse (opposite NES sign) or mimic (same sign), both padj < alpha.
gsea_concordance <- function(g, alpha = 0.05) {
  bind_rows(lapply(c("aDR", "YMP"), \(n) {
    j <- inner_join(g$aging %>% select(pathway, NES_age = NES, padj_age = padj), g[[n]] %>% select(pathway, NES, padj), by = "pathway")
    s <- filter(j, padj_age < alpha)
    tibble(intervention = n, aging_sets = nrow(s), significant_in_both = sum(s$padj < alpha),
           reversed = sum(s$padj < alpha & sign(s$NES) != sign(s$NES_age)), same_direction = sum(s$padj < alpha & sign(s$NES) == sign(s$NES_age)),
           NES_spearman_all_sets = round(cor(j$NES_age, j$NES, method = "spearman"), 2))
  }))
}

#' log2FC / padj of named genes in the aging and intervention contrasts.
gene_table <- function(res, symbols) {
  pick <- \(r, tag) r %>% filter(symbol %in% symbols) %>% transmute(symbol, !!paste0("lfc_", tag) := round(log2FoldChange, 2), !!paste0("padj_", tag) := signif(padj, 2))
  Reduce(\(a, b) left_join(a, b, by = "symbol"), list(pick(res$aging, "aging"), pick(res$aDR, "aDR"), pick(res$YMP, "YMP")))
}

# ---- team common outputs (results/common/, same file names and columns for every teammate) ----------------

#' Common reversal measure (team setup). Genes = aging DEGs (padj < alpha) tested in the treatment comparison.
#' spearman_rho: aging vs treatment log2FC (shrunken), negative = reversal; pct_reversed: % with opposite sign, binomial test vs 50 %.
reversal_measure <- function(aging, interv, alpha = 0.05) {
  a <- filter(aging, !is.na(padj), padj < alpha)
  j <- tibble(age = lfc_of(a), int = lfc_of(interv)[match(a$gene_id, interv$gene_id)]) %>% filter(!is.na(int))
  ct <- suppressWarnings(cor.test(j$age, j$int, method = "spearman"))
  rev <- sum(sign(j$age) != sign(j$int))
  tibble(n_genes = nrow(j), spearman_rho = round(unname(ct$estimate), 3), spearman_p = signif(ct$p.value, 2),
         pct_reversed = round(100 * rev / nrow(j), 1), binom_p = signif(binom.test(rev, nrow(j), 0.5)$p.value, 2))
}

#' Write the team's common files: de_<comparison>.csv, gsea_hallmark.csv, gsea_gobp.csv, reversal.csv.
#' res = list(aging, aDR, YMP) of DE tables with lfc_shrunk; gs_h / gs_go = lists of run_gsea() results per comparison.
write_common <- function(res, gs_h, gs_go, out = file.path(ROOT, "results/common")) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  iwalk(res, \(r, n) write_csv(select(r, gene_id, symbol, baseMean, log2FoldChange, lfc_shrunk, stat, pvalue, padj), file.path(out, paste0("de_", n, ".csv"))))
  g <- \(x) bind_rows(x, .id = "comparison") %>% select(comparison, pathway, NES, pval, padj, size)
  write_csv(g(gs_h), file.path(out, "gsea_hallmark.csv")); write_csv(g(gs_go), file.path(out, "gsea_gobp.csv"))
  rv <- bind_rows(lapply(c("aDR", "YMP"), \(n) reversal_measure(res$aging, res[[n]]) %>% mutate(treatment = n, .before = 1)))
  write_csv(rv, file.path(out, "reversal.csv"))
  rv
}
