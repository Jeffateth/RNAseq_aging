# TA check: does trimming change STAR's uniquely mapped proportion? Reads the two arms of
# results/trim_check/<SRR>/ (raw, trimmed; written by euler/06_trim_check.sh) and writes
# results/compare/trim_check.csv plus results/compare/trim_check.png.
suppressPackageStartupMessages(library(ggplot2))
root <- "results/trim_check"
srr  <- basename(list.dirs(root, recursive = FALSE))[1]
star_log <- function(arm) {
  f <- file.path(root, srr, arm, "Log.final.out")
  v <- readLines(f); get <- function(k) as.numeric(gsub("[%\t ]", "", sub(".*\\|\t", "", v[grep(k, v, fixed = TRUE)][1])))
  data.frame(arm = arm, input_pairs = get("Number of input reads"),
             mean_input_len = get("Average input read length"), unique_pct = get("Uniquely mapped reads %"),
             multi_pct = get("% of reads mapped to multiple loci"),
             too_short_pct = get("% of reads unmapped: too short"), mean_mapped_len = get("Average mapped length"))
}
tab <- do.call(rbind, lapply(c("raw", "trimmed"), star_log))
tab$arm <- factor(tab$arm, levels = c("raw", "trimmed"),
                  labels = c("no fastp", "pipeline fastp"))
dir.create("results/compare", showWarnings = FALSE, recursive = TRUE)
write.csv(tab, "results/compare/trim_check.csv", row.names = FALSE)
print(tab)

p <- ggplot(tab, aes(arm, unique_pct)) +
  geom_col(width = 0.5, fill = "grey35") +
  geom_text(aes(label = sprintf("%.2f %%", unique_pct)), vjust = -0.4, size = 5) +
  coord_cartesian(ylim = c(min(tab$unique_pct) - 2, max(tab$unique_pct) + 2)) +
  labs(title = paste0("STAR uniquely mapped reads, ", srr),
       subtitle = "same index and STAR call, only the reads differ",
       x = NULL, y = "uniquely mapped (%)") +
  theme_bw(base_size = 16)
ggsave("results/compare/trim_check.png", p, width = 7, height = 5, dpi = 150)
