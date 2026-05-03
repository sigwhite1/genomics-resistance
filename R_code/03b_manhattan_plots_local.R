# =============================================================================
# 03b_manhattan_plots_local.R
# Purpose: Local plot tweaking using pre-computed RDS files from Savio
# Run this in RStudio after transferring RDS files from Savio
# =============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)
library(gridExtra)
library(grid)

# Load pre-computed data from Savio
manhattan_plot_data <- readRDS("~/Desktop/IL Genomics/manhattan_plot_data.rds")
pca_plot_data       <- readRDS("~/Desktop/IL Genomics/pca_plot_data.rds")
qq_wza              <- readRDS("~/Desktop/IL Genomics/qq_wza.rds")
qq_pca_data         <- readRDS("~/Desktop/IL Genomics/qq_pca_data.rds")
top10               <- readRDS("~/Desktop/IL Genomics/top10.rds")

out_dir <- "~/Desktop/IL Genomics/figures"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# FDR thresholds — hardcoded from Savio run
fdr_threshold_wza  <- -log10(0.05 / 571)   # 571 WZA windows
fdr_threshold_logP <- -log10(0.05 / 571)   # 571 PCA windows

# =============================================================================
# WZA Manhattan plot
# =============================================================================

p_wza <- ggplot(manhattan_plot_data, aes(x=position, y=logP, color=color)) +
  geom_point(alpha=0.75, size=1.5, na.rm=TRUE) +
  geom_hline(yintercept=fdr_threshold_wza, linetype="dashed",
             color="red", linewidth=1) +
  facet_grid(~chrom, scales="fixed", space="fixed") +
  scale_color_manual(values=c("0"="gray30","1"="gray60","padding"=NA),
                     guide="none") +
  labs(title="WZA scores by chromosome",
       x="Chromosome", y="-log10(WZA P-value)") +
  theme_minimal() +
  theme(
    strip.text=element_text(size=5, angle=90, hjust=1, vjust=0.5),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.spacing=unit(0.1, "lines")
  )

# =============================================================================
# PCA-LD50 Manhattan plot
# =============================================================================

p_pca <- ggplot(pca_plot_data, aes(x=position, y=logP, color=color)) +
  geom_point(alpha=0.75, size=1.5, na.rm=TRUE) +
  geom_hline(yintercept=fdr_threshold_logP, linetype="dashed",
             color="red", linewidth=1) +
  facet_grid(~chrom, scales="fixed", space="fixed") +
  scale_color_manual(values=c("0"="gray30","1"="gray60","padding"=NA),
                     guide="none") +
  labs(title="Manhattan plot of PCA-LD50 correlation",
       x="Chromosome", y=expression(-log[10](italic(P)))) +
  theme_minimal() +
  theme(
    strip.text=element_text(size=5, angle=90, hjust=1, vjust=0.5),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.spacing=unit(0.1, "lines")
  )

# =============================================================================
# QQ plots
# =============================================================================

p_qq_wza <- ggplot(qq_wza, aes(x=expected, y=observed)) +
  geom_point(size=1.5, alpha=0.7, color="blue") +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red") +
  theme_minimal() +
  labs(title="WZA QQ plot",
       x="Expected -log10(P)", y="Observed -log10(P)")

p_qq_pca <- ggplot(qq_pca_data, aes(x=Expected, y=Observed)) +
  geom_point(alpha=0.7, size=2, color="blue") +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
  theme_minimal() +
  labs(title="PCA-LD50 QQ plot",
       x="Expected -log10(P)", y="Observed -log10(P)")

# =============================================================================
# Save figures
# =============================================================================

out_pdf <- file.path(out_dir, "Figure_5_manhattan_wza.pdf")
pdf(out_pdf, width=14, height=10)

# Page 1: Combined Manhattan
print(p_pca / p_wza + plot_layout(heights=c(1,1)))

# Page 2: QQ plots
print(p_qq_pca + p_qq_wza)

# Page 3: Top 10 table
grid.newpage()
grid.draw(tableGrob(top10, rows=NULL))
grid.text("Top 10 Genomic Windows by WZA",
          y=unit(1,"npc") - unit(2,"mm"),
          gp=gpar(fontsize=14, fontface="bold"))

dev.off()
cat("Figures saved to", out_pdf, "\n")
