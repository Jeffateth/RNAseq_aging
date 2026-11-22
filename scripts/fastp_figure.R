# fastp before / after trimming for one run (read 1): base content and quality per position, large text for slides.
# Same curves as the fastp HTML report, redrawn from results/qc/<SRR>.fastp.json.
# Run from the repository root: .pixi/envs/default/bin/Rscript scripts/fastp_figure.R [SRR]   (default SRR21355385)
# Writes results/qc/fastp_before_after_<SRR>.png
suppressPackageStartupMessages({library(jsonlite); library(tidyverse)})
args <- commandArgs(trailingOnly = TRUE)
srr <- if (length(args)) args[1] else "SRR21355385"
j <- fromJSON(file.path("results/qc", paste0(srr, ".fastp.json")))

curves <- function(stage, what) {
  x <- j[[paste0("read1_", stage, "_filtering")]][[what]]
  bind_rows(lapply(names(x), \(b) tibble(base = b, position = seq_along(x[[b]]), value = unlist(x[[b]]))))
}
d <- bind_rows(
  curves("before", "content_curves") %>% mutate(stage = "Before trimming", panel = "Base content (%)", value = 100 * value),
  curves("after",  "content_curves") %>% mutate(stage = "After trimming",  panel = "Base content (%)", value = 100 * value),
  curves("before", "quality_curves") %>% mutate(stage = "Before trimming", panel = "Quality (Phred)"),
  curves("after",  "quality_curves") %>% mutate(stage = "After trimming",  panel = "Quality (Phred)")) %>%
  filter(base %in% c("A", "T", "C", "G", "GC", "mean")) %>%
  mutate(base = factor(ifelse(base %in% c("GC", "mean"), "black", base), c("A", "T", "C", "G", "black")),
         stage = factor(stage, c("Before trimming", "After trimming")))

cols <- c(A = "#1baf7a", T = "#a35bd6", C = "#eb6834", G = "#2a78d6", black = "#222222")
# poly-G: G rises towards the 3' end before trimming only
annotate_poly <- geom_text(data = tibble(stage = factor("Before trimming", levels(d$stage)), panel = "Base content (%)", position = 112, value = 36),
                           aes(position, value, label = "poly-G"), inherit.aes = FALSE, size = 13, colour = "#2a78d6", fontface = "bold")
p <- ggplot(d, aes(position, value, colour = base)) +
  geom_line(linewidth = 1.6) +
  facet_grid(panel ~ stage, scales = "free_y", switch = "y") +
  scale_colour_manual(values = cols, labels = c(A = "A", T = "T", C = "C", G = "G", black = "GC (top) / mean (bottom)")) +
  annotate_poly +
  scale_x_continuous(breaks = c(1, 50, 100, 150)) +
  labs(x = "Position in read 1 (bases)", y = NULL, colour = NULL) +
  theme_minimal(base_size = 40) +
  theme(panel.grid.minor = element_blank(), strip.placement = "outside", strip.text = element_text(size = rel(0.95)),
        legend.position = "bottom", legend.key.width = unit(3, "lines"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0.3, "lines"), panel.spacing = unit(2, "lines"),
        plot.margin = margin(4, 10, 4, 4)) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 4)))
ggsave(file.path("results/qc", paste0("fastp_before_after_", srr, ".png")), p, width = 22, height = 12, dpi = 150, bg = "white")

s <- j$summary; f <- j$filtering_result; n <- s$before_filtering$total_reads
cat(sprintf("%s: %.1f M reads; kept %.1f %%; low quality %.1f %%; too short %.2f %%; Q30 %.1f -> %.1f %%; GC %.1f -> %.1f %%\n",
            srr, n / 1e6, 100 * f$passed_filter_reads / n, 100 * f$low_quality_reads / n, 100 * f$too_short_reads / n,
            100 * s$before_filtering$q30_rate, 100 * s$after_filtering$q30_rate, 100 * s$before_filtering$gc_content, 100 * s$after_filtering$gc_content))
