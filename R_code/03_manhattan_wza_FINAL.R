# =============================================================================
# 03_manhattan_wza_FINAL.R
#
# Purpose: Generate Manhattan plots, QQ plots, and the convergent-window
#          table for Figure 4, reading directly from the single verified
#          source of truth (joint_results_FINAL.csv / convergent_windows_FINAL.csv,
#          produced by wza_final.R) rather than recomputing WZA/PCA from
#          wza_cache/ and snprelate_1mb.csv a second time.
#
# WHY THIS SCRIPT EXISTS:
#   03_manhattan_wza_CORRECTED.R correctly fixed the sample-mismatch bug,
#   the WZA rank-tie bug, and the WZA Manhattan y-axis bug -- but it
#   recomputes WZA normalization independently using the OLD
#   quadratic-regression dispersion model (lm(ZW ~ poly(n_SNPs,2))) rather
#   than the Gamma GLM fix, and applies no window-inclusion filter
#   (full-size + n_SNPs >= 100). That produced a figures.pdf containing a
#   spurious |Z_W| = 161.72 outlier from a 90-SNP window that should have
#   been excluded, and a convergent-windows table built on the wrong
#   denominator. wza_final.R has the complete, correct pipeline and writes
#   joint_results_FINAL.csv / convergent_windows_FINAL.csv -- this script
#   only handles presentation (plots + table) of that already-verified data,
#   with no independent recomputation of anything.
#
# Prerequisites (all on Savio, both written by wza_final.R):
#   - joint_results_FINAL.csv
#   - convergent_windows_FINAL.csv
#
# Outputs:
#   - figures.pdf  (Manhattan plots, QQ plots, convergent-window table)
# =============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)
library(gridExtra)
library(grid)

analysis_dir <- "/global/scratch/users/pswhite/analysis"
out_pdf      <- file.path(analysis_dir, "figures.pdf")

joint <- read.csv(file.path(analysis_dir, "joint_results_FINAL.csv"))
conv  <- read.csv(file.path(analysis_dir, "convergent_windows_FINAL.csv"))

cat("Loaded joint_results_FINAL.csv:", nrow(joint), "windows\n")
cat("Max |Z_W,norm|:", round(max(joint$abs_ZW, na.rm = TRUE), 3), "\n")
cat("Windows with n_SNPs < 100:", sum(joint$n_SNPs < 100, na.rm = TRUE), "\n")
cat("Convergent windows:", nrow(conv), "\n")

# =============================================================================
# 1. PARSE CHROMOSOME / POSITION FROM clean_window
# =============================================================================
# NOTE: this assumes clean_window looks like "NC_071295.1_6750000_7000000.vcf"
# (no "window_" / "slidingwindow_window_" prefix -- both were already stripped
# upstream in wza_final.R). The printed sanity check below confirms this
# before anything is plotted; if chrom/position come out NA, the format
# differs from what's assumed here and the regex needs adjusting.

parse_window <- function(df) {
  df %>%
    mutate(
      chrom        = sub("((NC|NW|CM)_[0-9]+\\.[0-9]+)_.*", "\\1", clean_window),
      position     = as.numeric(sub(".*_([0-9]+)_([0-9]+)\\.vcf$", "\\1", clean_window)),
      end_position = as.numeric(sub(".*_([0-9]+)_([0-9]+)\\.vcf$", "\\2", clean_window))
    ) %>%
    mutate(
      chrom_index = as.numeric(factor(chrom, levels = unique(chrom))),
      color       = as.factor(chrom_index %% 2)
    )
}

manhattan_data <- parse_window(joint)

cat("\n--- SANITY CHECK: confirm parsing worked before trusting the plots ---\n")
cat("Distinct chromosomes found:", length(unique(manhattan_data$chrom)), "(expect 31)\n")
cat("Rows with NA position after parsing:", sum(is.na(manhattan_data$position)), "(expect 0)\n")
print(head(manhattan_data[, c("clean_window", "chrom", "position", "end_position")]))
cat("--- end sanity check ---\n\n")

# =============================================================================
# 2. WZA MANHATTAN (y = |Z_W,norm|)
# =============================================================================

wza_z_threshold <- qnorm(1 - (0.05 / nrow(manhattan_data)) / 2)

global_max_pos <- max(manhattan_data$end_position, na.rm = TRUE)
padding_wza <- manhattan_data %>%
  distinct(chrom, chrom_index) %>%
  mutate(position = global_max_pos, end_position = global_max_pos,
         abs_ZW = NA, color = "padding", clean_window = NA)

manhattan_plot_data <- bind_rows(manhattan_data, padding_wza) %>%
  mutate(color = factor(ifelse(is.na(as.character(color)), "padding", as.character(color)),
                         levels = c("0", "1", "padding")))

p_wza <- ggplot(manhattan_plot_data, aes(x = position, y = abs_ZW, color = color)) +
  geom_point(alpha = 0.75, size = 1.2, na.rm = TRUE) +
  geom_hline(yintercept = wza_z_threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
  facet_grid(~chrom, scales = "fixed", space = "fixed") +
  scale_color_manual(values = c("0" = "gray30", "1" = "gray60", "padding" = NA), guide = "none") +
  labs(title = "WZA scores by chromosome",
       x = "Chromosome", y = expression(group("|", Z[W] ~ norm, "|"))) +
  theme_minimal() +
  theme(strip.text = element_text(size = 5, angle = 45, hjust = 1, vjust = 1),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

# =============================================================================
# 3. PCA-LD50 MANHATTAN (y = -log10(p))
# =============================================================================

manhattan_data_pca <- manhattan_data %>% mutate(logP = -log10(p_value))
bonferroni_threshold_logP <- -log10(0.05 / nrow(manhattan_data_pca))

padding_pca <- manhattan_data_pca %>%
  distinct(chrom, chrom_index) %>%
  mutate(position = global_max_pos, end_position = global_max_pos,
         logP = NA, color = "padding", clean_window = NA)

pca_plot_data <- bind_rows(manhattan_data_pca, padding_pca) %>%
  mutate(color = factor(ifelse(is.na(as.character(color)), "padding", as.character(color)),
                         levels = c("0", "1", "padding")))

p_pca <- ggplot(pca_plot_data, aes(x = position, y = logP, color = color)) +
  geom_point(alpha = 0.75, size = 1.2, na.rm = TRUE) +
  geom_hline(yintercept = bonferroni_threshold_logP, linetype = "dashed", color = "red", linewidth = 0.8) +
  facet_grid(~chrom, scales = "fixed", space = "fixed") +
  scale_color_manual(values = c("0" = "gray30", "1" = "gray60", "padding" = NA), guide = "none") +
  labs(title = "Manhattan plot of PCA-LD50 correlation",
       x = "Chromosome", y = expression(-log[10](italic(P)))) +
  theme_minimal() +
  theme(strip.text = element_text(size = 5, angle = 45, hjust = 1, vjust = 1),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

# =============================================================================
# 4. QQ PLOTS
# =============================================================================

qq_wza <- joint %>%
  filter(!is.na(abs_ZW)) %>%
  arrange(desc(abs_ZW)) %>%
  mutate(observed = abs_ZW, expected = qnorm(1 - (ppoints(n()) / 2)))

p_qq_wza <- ggplot(qq_wza, aes(x = expected, y = observed)) +
  geom_point(size = 1.2, alpha = 0.7, color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "WZA QQ plot",
       x = expression(Expected ~ group("|", Z, "|")),
       y = expression(Observed ~ group("|", Z, "|")))

qq_pca_data <- data.frame(
  Expected = -log10(ppoints(nrow(joint))),
  Observed = -log10(sort(joint$p_value))
)

p_qq_pca <- ggplot(qq_pca_data, aes(x = Expected, y = Observed)) +
  geom_point(alpha = 0.7, size = 1.5, color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "PCA-LD50 QQ plot",
       x = expression(Expected ~ -log[10](italic(P))),
       y = expression(Observed ~ -log[10](italic(P))))

# =============================================================================
# 5. CONVERGENT WINDOWS TABLE
# =============================================================================

conv <- conv[order(conv$pve_rank), ]

conv_display <- data.frame(
  Window     = conv$clean_window,
  Tier       = conv$tier,
  r          = round(conv$correlation, 3),
  PVE        = round(conv$PVE, 3),
  `PVE rank` = conv$pve_rank,
  Z_W        = round(conv$ZW_norm, 2),
  `WZA rank` = conv$wza_rank,
  check.names = FALSE
)

cat("\nConvergent windows (from convergent_windows_FINAL.csv):\n")
print(conv_display, row.names = FALSE)

# =============================================================================
# 6. SAVE TO PDF
# =============================================================================

cat("\nSaving figures to", out_pdf, "...\n")
pdf(out_pdf, width = 14, height = 10)

print(p_pca / p_wza + plot_layout(heights = c(1, 1)))
print(p_qq_pca + p_qq_wza)

grid.newpage()
grid.draw(tableGrob(conv_display, rows = NULL, theme = ttheme_default(base_size = 8)))
grid.text("Convergent genomic windows (top 1% and top 2% by both PVE and WZA)",
          y = unit(1, "npc") - unit(2, "mm"),
          gp = gpar(fontsize = 13, fontface = "bold"))

dev.off()
cat("Done. Figures saved to", out_pdf, "\n")
cat("(Old figures.pdf has been overwritten. If you want to compare against the\n")
cat(" stale version first, rename or copy it before running this script.)\n")
