# Why DESeq2 instead of a per-gene ANOVA/ANCOVA (~ sex + age)? One figure, aging analysis (old vs young, ~ sex + age, 11 mice).
# Run from the repository root (seconds): .pixi/envs/default/bin/Rscript scripts/deseq2_vs_ancova.R
# Writes results/diagnostics/deseq2_vs_ancova.{png,pdf}; DEG tables are unchanged.
setwd("notebooks"); suppressMessages(source("../R/helpers.R")); 
out <- file.path(ROOT, "results/diagnostics"); dir.create(out, showWarnings = FALSE)

samples <- load_samples()
s_age <- samples %>% filter(cohort == "aging", rownames(.) != "SRR21355389")
d <- make_dds(load_rsem(rownames(s_age)), s_age, ~ sex + age_class, min_n = 4)   # same model as notebook 03
cd <- as.data.frame(colData(d))
y <- log2(counts(d, normalized = TRUE) + 1)

# ANOVA/ANCOVA: one linear model per gene on log2 normalized counts, same design (~ sex + age), t-test, BH
ancova <- function(age) {
  X <- model.matrix(~ cd$sex + age); f <- lm.fit(X, t(y)); df <- nrow(X) - ncol(X)
  sd <- sqrt(colSums(f$residuals^2) / df)
  list(padj = p.adjust(2 * pt(-abs(f$coefficients[3, ] / (sd * sqrt(solve(crossprod(X))[3, 3]))), df), "BH"), sd = sd)
}

# Real labels
a <- ancova(cd$age_class); r <- results(d, contrast = c("age_class", "old", "young"))
deg_d <- which(r$padj < 0.05); deg_a <- which(a$padj < 0.05)

# One plot: aging DEGs per method, split into genes found by both and by one method only
shared <- length(intersect(deg_a, deg_d))
pd <- tibble(method = factor(c("ANCOVA\n(two-way ANOVA)", "ANCOVA\n(two-way ANOVA)", "DESeq2", "DESeq2"), c("DESeq2", "ANCOVA\n(two-way ANOVA)")),
             part = factor(c("found by both", "only this method", "found by both", "only this method"),
                           c("only this method", "found by both")),
             n = c(shared, length(deg_a) - shared, shared, length(deg_d) - shared))
tot <- tibble(method = factor(c("ANCOVA\n(two-way ANOVA)", "DESeq2"), c("DESeq2", "ANCOVA\n(two-way ANOVA)")), n = c(length(deg_a), length(deg_d)))
fig <- ggplot(pd, aes(n, method, fill = part)) +
  geom_col(width = 0.6) +
  geom_text(data = tot, aes(n, method, label = n), inherit.aes = FALSE, hjust = -0.25, size = 18, fontface = "bold") +
  annotate("text", x = shared / 2, y = c(1, 2), label = paste(shared, "shared"), size = 8, colour = "white", fontface = "bold") +
  scale_fill_manual(values = c("found by both" = "#1baf7a", "only this method" = "#2a78d6"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = NULL, title = "Aging genes found (old vs young)",
       subtitle = "Same 11 mice, same design (age + sex), padj < 0.05") +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 32),
        plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 26, colour = "grey30"))
ggsave(file.path(out, "deseq2_vs_ancova.png"), fig, width = 16, height = 7, dpi = 100, bg = "white")
ggsave(file.path(out, "deseq2_vs_ancova.pdf"), fig, width = 16, height = 7, bg = "white")
cat(sprintf("DESeq2 %d, ANCOVA %d, shared %d\n", length(deg_d), length(deg_a), shared))
