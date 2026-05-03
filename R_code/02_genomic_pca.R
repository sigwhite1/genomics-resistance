# =============================================================================
# 02_genomic_pca.R
#
# Purpose: Genomic PCA analyses and figures
#   - Figure S2: Neighbor-joining tree (tips colored by LD50)
#   - Figure S: Genomic PCA with stock population (two stock clusters)
#   - Figure 3: Genomic PC1 vs LD50 (headline result, R = -0.64, p = 0.026)
#
# Inputs:
#   - ~/Desktop/Trade-offs Data/Updated Analyses/genetic_distances.csv
#   - ~/Desktop/Trade-offs Data/stock_pca_results.csv
#   - ~/Desktop/Trade-offs Data/ld50_calculated.csv
#
# Outputs:
#   - Figure_S1_NJ_tree.pdf
#   - Figure_S2_stock_PCA.pdf
#   - Figure_3_genomic_PC1_vs_LD50.pdf
# =============================================================================

library(ape)
library(phangorn)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)

# ── USER SETTINGS ─────────────────────────────────────────────────────────────
out_dir <- "~/Desktop/IL Genomics/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# ──────────────────────────────────────────────────────────────────────────────

# ── Shared data ───────────────────────────────────────────────────────────────

# LD50 values calculated in 01_life_history.R via MASS::dose.p
# (do NOT hardcode here — run 01_life_history.R first to generate this file)
ld50_file <- "~/Desktop/ld50_calculated.csv"
if (!file.exists(ld50_file)) stop("ld50_calculated.csv not found. Run 01_life_history.R first.")
ld50_table  <- read.csv(ld50_file)
ld50_values <- setNames(ld50_table$Avg_LD50, ld50_table$Population)
# Exclude Stock from inbred-line analyses
ld50_values <- ld50_values[names(ld50_values) != "Stock"]

# Color ramp for LD50
color_ramp <- colorRampPalette(c("steelblue", "gold", "firebrick"))
n_colors   <- 100
colors     <- color_ramp(n_colors)

# BAM path → IL name mapping
il_map <- c(
  "IL-1"="IL-1",  "IL-2"="IL-2",  "IL-3"="IL-4",  "IL-4"="IL-7",
  "IL-5"="IL-9",  "IL-6"="IL-10", "IL-7"="IL-14", "IL-8"="IL-17",
  "IL-9"="IL-19", "IL-10"="IL-24","IL-11"="IL-29","IL-12"="IL-36"
)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

# Distance matrix (for NJ tree)
dist_mat      <- read.csv("~/Desktop/genetic_distances.csv", row.names=1)
dist_mat_full <- as.matrix(dist_mat)
dist_obj      <- as.dist(dist_mat_full)

# PCA results (inbred + stock)
pca_raw <- read.csv("~/Desktop/stock_pca_results.csv")

# Clean sample names
pca <- pca_raw %>%
  mutate(sample_clean = case_when(
    grepl("I10\\.bam", sample) ~ "IL-10",
    grepl("I11\\.bam", sample) ~ "IL-11",
    grepl("I12\\.bam", sample) ~ "IL-12",
    grepl("I1\\.bam",  sample) ~ "IL-1",
    grepl("I2\\.bam",  sample) ~ "IL-2",
    grepl("I3\\.bam",  sample) ~ "IL-3",
    grepl("I4\\.bam",  sample) ~ "IL-4",
    grepl("I5\\.bam",  sample) ~ "IL-5",
    grepl("I6\\.bam",  sample) ~ "IL-6",
    grepl("I7\\.bam",  sample) ~ "IL-7",
    grepl("I8\\.bam",  sample) ~ "IL-8",
    grepl("I9\\.bam",  sample) ~ "IL-9",
    grepl("S[0-9]",   sample)  ~ sub(".*(S[0-9]+).*", "\\1", sample),
    TRUE ~ sample
  )) %>%
  mutate(
    label = ifelse(population == "Inbred", il_map[sample_clean], sample_clean),
    LD50  = ifelse(population == "Inbred", ld50_values[label], NA)
  )

# =============================================================================
# 2. FIGURE S2 — Neighbor-joining tree (all 12 inbred lines)
# =============================================================================

tree <- nj(dist_obj)

# Use same gradient as other figures
color_ramp_fig <- colorRampPalette(c("steelblue", "gold", "firebrick"))
n_colors       <- 100
colors_fig     <- color_ramp_fig(n_colors)

tip_ld50    <- ld50_values[tree$tip.label]
ld50_scaled <- as.integer((tip_ld50 - min(tip_ld50)) /
                            (max(tip_ld50) - min(tip_ld50)) * (n_colors - 1)) + 1
tip_colors  <- colors_fig[ld50_scaled]

pdf("~/Desktop/figures/Figure_S2_NJ_tree.pdf", width=7, height=8)
par(mar=c(5, 1, 2, 5))

plot(tree,
     type="phylogram", tip.color=tip_colors,
     edge.width=2.5, cex=1.2, font=3,
     direction="rightwards", no.margin=FALSE,
     x.lim=c(0, 420))

add.scale.bar(x=0, y=0.8, cex=0.85, col="gray40")
text(x=55, y=0.3, "50 allelic differences", cex=0.8, col="gray40")

# Italicized title
title(expression(atop(
  "Neighbor-joining tree of genomic relatedness",
  paste("among inbred lines of ", italic("Plodia interpunctella"))
)), cex.main=1.0, font.main=1)

# Gradient legend — moved down and condensed
legend_x      <- 290
legend_y_top  <- 7     # lowered from 12
legend_y_bot  <- 3     # lowered from 7, condensed height
legend_width  <- 15
n_legend      <- 100

legend_colors <- color_ramp_fig(n_legend)
y_steps       <- seq(legend_y_bot, legend_y_top, length.out=n_legend + 1)

for (i in 1:n_legend) {
  rect(xleft   = legend_x,
       xright  = legend_x + legend_width,
       ybottom = y_steps[i],
       ytop    = y_steps[i + 1],
       col     = legend_colors[i],
       border  = NA)
}

rect(xleft   = legend_x,
     xright  = legend_x + legend_width,
     ybottom = legend_y_bot,
     ytop    = legend_y_top,
     col     = NA, border="grey30", lwd=0.8)

ld50_min <- round(min(tip_ld50), 2)
ld50_max <- round(max(tip_ld50), 2)
ld50_mid <- round(mean(c(ld50_min, ld50_max)), 2)

text(x=legend_x + legend_width + 2, y=legend_y_bot,
     labels=ld50_min, adj=0, cex=0.8)
text(x=legend_x + legend_width + 2, y=mean(c(legend_y_bot, legend_y_top)),
     labels=ld50_mid, adj=0, cex=0.8)
text(x=legend_x + legend_width + 2, y=legend_y_top,
     labels=ld50_max, adj=0, cex=0.8)

text(x=legend_x + legend_width / 2, y=legend_y_top + 0.5,
     labels="LD50", adj=0.5, cex=0.9, font=2)

dev.off()

# =============================================================================
# EXPLORATORY — Inbred-only PCA (PC1 vs PC2, coloured by LD50)
# =============================================================================

inbred_pca <- pca %>%
  filter(population == "Inbred") %>%
  mutate(LD50 = ld50_values[label])

fig_inbred_pca <- ggplot(inbred_pca, aes(x=PC1, y=PC2)) +
  geom_point(aes(fill=LD50), size=5, shape=21, stroke=0.5) +
  scale_fill_gradientn(colors=c("steelblue","gold","firebrick"),
                       name="LD50") +
  geom_text_repel(aes(label=label), fontface="italic",
                  size=3.5, box.padding=0.4, seed=42) +
  theme_classic() +
  labs(title="Genomic PCA: inbred lines only",
       x=paste0("PC1 (12.8%)"),
       y=paste0("PC2 (4.8%)")) +
  theme(plot.title=element_text(size=12))

ggsave(file.path(out_dir, "FigureS_inbred_PCA.pdf"),
       fig_inbred_pca, width=7, height=6)
cat("Inbred-only PCA saved.\n")


# =============================================================================
# 3. FIGURE 2 — Genomic PCA with stock population (two-panel)
# Panel A: Full PCA showing inbred lines and both stock clusters
# Panel B: Zoomed inbred-only view with clearer labels
# =============================================================================

# Lock cluster colors explicitly so they don't swap between runs
# Check which cluster has higher mean PC2 (cluster 1 = upper = orange)
cluster_means <- pca %>%
  filter(population == "Stock") %>%
  group_by(cluster) %>%
  summarise(mean_PC2 = mean(PC2), .groups = "drop") %>%
  arrange(desc(mean_PC2))

upper_cluster <- cluster_means$cluster[1]  # higher PC2 = cluster 1 = orange
lower_cluster <- cluster_means$cluster[2]  # lower PC2 = cluster 2 = purple

# Panel A: Full plot
panel_a <- ggplot(pca, aes(x=PC1, y=PC2)) +
  geom_point(data=pca[pca$cluster == upper_cluster, ],
             color="darkorange", size=1.5, alpha=0.7) +
  geom_point(data=pca[pca$cluster == lower_cluster, ],
             color="darkorchid", size=1.5, alpha=0.7) +
  geom_point(data=pca[pca$cluster == "Inbred line", ],
             aes(fill=LD50), size=4, shape=21, stroke=0.5) +
  scale_fill_gradientn(colors=c("steelblue","gold","firebrick"),
                       name="LD50") +
  annotate("text", x=50, y=210, label="Stock cluster 1",
           color="darkorange", size=3.5, fontface="bold") +
  annotate("text", x=50, y=-210, label="Stock cluster 2",
           color="darkorchid", size=3.5, fontface="bold") +
  annotate("rect",
           xmin=min(pca$PC1[pca$cluster == "Inbred line"]) - 20,
           xmax=max(pca$PC1[pca$cluster == "Inbred line"]) + 20,
           ymin=min(pca$PC2[pca$cluster == "Inbred line"]) - 20,
           ymax=max(pca$PC2[pca$cluster == "Inbred line"]) + 20,
           color="grey40", fill=NA, linetype="dashed", linewidth=0.4) +
  theme_classic() +
  labs(tag="A",
       x="PC1 (12.8%)", y="PC2 (4.8%)") +
  theme(legend.position="none")

# Panel B: Zoomed inbred-only view
inbred_only <- pca %>%
  filter(population == "Inbred") %>%
  mutate(LD50 = ld50_values[label])

# Compute axis limits with tighter padding
x_pad <- diff(range(inbred_only$PC1)) * 0.08
y_pad <- diff(range(inbred_only$PC2)) * 0.15

panel_b <- ggplot(inbred_only, aes(x=PC1, y=PC2)) +
  geom_point(aes(fill=LD50), size=4, shape=21, stroke=0.5) +
  scale_fill_gradientn(colors=c("steelblue","gold","firebrick"),
                       name="LD50") +
  geom_text_repel(aes(label=label), fontface="italic",
                  size=3.5, box.padding=0.5,
                  max.overlaps=20, seed=42) +
  coord_cartesian(
    xlim=c(min(inbred_only$PC1) - x_pad, max(inbred_only$PC1) + x_pad),
    ylim=c(min(inbred_only$PC2) - y_pad, max(inbred_only$PC2) + y_pad)
  ) +
  theme_classic() +
  labs(tag="B",
       x="PC1 (12.8%)", y="PC2 (4.8%)") +
  theme(legend.position="right")

# Combine panels with shared legend
fig2_combined <- panel_a + panel_b +
  plot_layout(widths=c(1, 1), guides="collect") &
  theme(legend.position="right")

ggsave(file.path(out_dir, "Figure_2_genomic_PCA_twopanel.pdf"),
       fig3_combined, width=14, height=6)
cat("Figure 2 (two-panel) saved.\n")

# =============================================================================
# 4. FIGURE 3 — Genomic PC1 vs LD50 (headline result)
# =============================================================================

inbred_only <- pca %>%
  filter(population == "Inbred") %>%
  mutate(LD50 = ld50_values[label])

cor_fig3 <- cor.test(inbred_only$PC1, inbred_only$LD50)
r_fig3   <- round(cor_fig3$estimate, 2)
p_fig3   <- round(cor_fig3$p.value, 3)

cat("\nGENOMIC PC1-LD50 CORRELATION:\n")
cat("r =", r_fig3, "  p =", p_fig3, "\n")

fig3 <- ggplot(inbred_only, aes(x=PC1, y=LD50)) +
  geom_smooth(method="lm", color="black", linetype="solid",
              se=TRUE, fill="grey85", linewidth=0.9) +
  geom_point(aes(fill=LD50), size=3.5, shape=21, stroke=0.5) +
  scale_fill_gradientn(colors=c("steelblue","gold","firebrick"),
                       guide="none") +
  geom_text_repel(aes(label=label), fontface="italic",
                  size=3, box.padding=0.4, seed=42,
                  nudge_y=ifelse(inbred_only$label == "IL-29", 0.025, 0)) +
  annotate("text", x=Inf, y=-Inf,
           label=paste0("r = ", r_fig3, ",  p = ", p_fig3),
           hjust=1.1, vjust=-0.5, size=4) +
  labs(title="Genomic structure predicts viral resistance",
       x="Genomic PC1 (12.8% variance explained)",
       y="LD50") +
  theme_classic()

ggsave(file.path(out_dir, "Figure_3_genomic_PC1_vs_LD50.pdf"), fig4,
       width=7, height=6)
cat("Figure 3 saved.\n")

# Sensitivity analysis: exclude IL-4
inbred_no_IL4 <- inbred_only %>% filter(label != "IL-4")
cor_no_IL4 <- cor.test(inbred_no_IL4$PC1, inbred_no_IL4$LD50)
r_no_IL4   <- round(cor_no_IL4$estimate, 2)
p_no_IL4   <- round(cor_no_IL4$p.value, 3)

cat("SENSITIVITY: Genomic PC1-LD50 excluding IL-4:\n")
cat("R =", r_no_IL4, "  p =", p_no_IL4, "\n")
