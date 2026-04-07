# =============================================================================
# 02_genomic_pca.R
#
# Purpose: Genomic PCA analyses and figures
#   - Figure S1: Neighbour-joining tree (tips colored by LD50)
#   - Figure S2: Genomic PCA with stock population (two stock clusters)
#   - Figure 3: Genomic PC1 vs LD50 (headline result, R = -0.64, p = 0.026)
#
# Inputs:
#   - ~/Desktop/Trade-offs Data/Updated Analyses (April 2026)/genetic_distances.csv
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

# ── USER SETTINGS ─────────────────────────────────────────────────────────────
out_dir <- "~/Desktop/IL Genomics/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# ──────────────────────────────────────────────────────────────────────────────

# ── Shared data ───────────────────────────────────────────────────────────────

# LD50 values calculated in 01_life_history.R via MASS::dose.p
# (do NOT hardcode here — run 01_life_history.R first to generate this file)
ld50_file <- "~/Desktop/IL Genomics/ld50_calculated.csv"
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
dist_mat      <- read.csv("~/Desktop/IL Genomics/Updated Analyses (April 2026)/genetic_distances.csv", row.names=1)
dist_mat_full <- as.matrix(dist_mat)
dist_obj      <- as.dist(dist_mat_full)

# PCA results (inbred + stock)
pca_raw <- read.csv("~/Desktop/IL Genomics/stock_pca_results.csv")

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
# 2. FIGURE S1 — Neighbor-joining tree (all 12 inbred lines)
# =============================================================================

tree <- nj(dist_obj)

tip_ld50    <- ld50_values[tree$tip.label]
ld50_scaled <- as.integer((tip_ld50 - min(tip_ld50)) /
                            (max(tip_ld50) - min(tip_ld50)) * (n_colors-1)) + 1
tip_colors  <- colors[ld50_scaled]

pdf(file.path(out_dir, "Figure_S1_NJ_tree.pdf"), width=7, height=8)
par(mar=c(5, 1, 2, 4))
plot(tree,
     type="phylogram", tip.color=tip_colors,
     edge.width=2.5, cex=1.2, font=3,
     direction="rightwards", no.margin=FALSE,
     x.lim=c(0, 420))
add.scale.bar(x=0, y=0.8, cex=0.85, col="gray40")
text(x=55, y=0.3, "50 allelic differences", cex=0.8, col="gray40")
legend(x=280, y=6,
       legend=c("Susceptible (low LD50)", "Intermediate", "Resistant (high LD50)"),
       fill=color_ramp(3), bty="n", cex=0.85)
title("Neighbor-joining tree of genomic relatedness\namong inbred lines of Plodia interpunctella",
      cex.main=1.0, font.main=1)
dev.off()
cat("Figure S1 saved.\n")

# =============================================================================
# 3. FIGURE S2 — Genomic PCA with stock population
# =============================================================================

# K-means clustering of stock individuals
stock_only    <- pca %>% filter(population == "Stock")
stock_clusters <- kmeans(stock_only[, c("PC1", "PC2")], centers=2, nstart=25)
stock_only$cluster <- as.factor(stock_clusters$cluster)

# Assign cluster labels to full dataset
pca <- pca %>%
  mutate(cluster = case_when(
    population == "Inbred" ~ "Inbred line",
    population == "Stock"  ~ paste0("Stock (cluster ",
                                    stock_clusters$cluster[
                                      match(sample_clean, stock_only$sample_clean)
                                    ], ")")
  ))

figS2 <- ggplot(pca, aes(x=PC1, y=PC2)) +
  geom_point(data=filter(pca, grepl("cluster 1", cluster)),
             color="darkorange", size=2, alpha=0.7) +
  geom_point(data=filter(pca, grepl("cluster 2", cluster)),
             color="darkorchid", size=2, alpha=0.7) +
  geom_point(data=filter(pca, cluster=="Inbred line"),
             aes(fill=LD50), size=4, shape=21, stroke=0.5) +
  scale_fill_gradientn(colors=c("steelblue","gold","firebrick"),
                       name="LD50") +
  geom_text_repel(data=filter(pca, cluster=="Inbred line"),
                  aes(label=label), size=3.5, fontface="italic",
                  max.overlaps=20, box.padding=0.4, seed=42) +
  annotate("text", x=50, y=210, label="Stock cluster 1",
           color="darkorange", size=3.5, fontface="bold") +
  annotate("text", x=50, y=-195, label="Stock cluster 2",
           color="darkorchid", size=3.5, fontface="bold") +
  theme_classic() +
  labs(title="Genomic PCA: inbred lines and stock population",
       subtitle="Stock individuals form two clusters consistent with two source populations",
       x="PC1 (12.8%)", y="PC2 (4.8%)") +
  theme(plot.title=element_text(size=12),
        plot.subtitle=element_text(size=10, color="gray40"))

ggsave(file.path(out_dir, "FigureS2_stock_PCA.pdf"), figS2,
       width=8, height=7)
cat("Figure S2 saved.\n")

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
# 4. FIGURE 3 — Genomic PC1 vs LD50 (headline result)
# =============================================================================

# Use inbred-only global PCA scores
inbred_only <- pca %>%
  filter(population == "Inbred") %>%
  mutate(LD50 = ld50_values[label])

cor_fig3 <- cor.test(inbred_only$PC1, inbred_only$LD50)
r_fig3   <- round(cor_fig3$estimate, 2)
p_fig3   <- round(cor_fig3$p.value, 3)

cat("\nGENOMIC PC1-LD50 CORRELATION:\n")
cat("R =", r_fig5, "  p =", p_fig5, "\n")
print(cor_fig5)

fig3 <- ggplot(inbred_only, aes(x=PC1, y=LD50)) +
  geom_smooth(method="lm", se=TRUE, color="gray70", fill="gray90") +
  geom_point(aes(fill=LD50), size=4, shape=21, stroke=0.5) +
  scale_fill_gradientn(colors=c("steelblue","gold","firebrick"),
                       guide="none") +
  geom_text_repel(aes(label=label), fontface="italic",
                  size=3.5, box.padding=0.4, seed=42) +
  annotate("text", x=Inf, y=-Inf,
           label=paste0("R = ", r_fig3, ",  p = ", p_fig3),
           hjust=1.1, vjust=-0.5, size=4, color="gray30") +
  theme_classic() +
  labs(x="Genomic PC1 (12.8% variance explained)",
       y="LD50",
       title="Genomic structure predicts viral resistance") +
  theme(plot.title=element_text(size=12))

ggsave(file.path(out_dir, "Figure_3_genomic_PC1_vs_LD50.pdf"), fig3,
       width=6, height=5)
cat("Figure 3 saved.\n")

# Sensitivity analysis: exclude IL-4
inbred_no_IL4 <- inbred_only %>% filter(label != "IL-4")

cor_no_IL4 <- cor.test(inbred_no_IL4$PC1, inbred_no_IL4$LD50)
r_no_IL4   <- round(cor_no_IL4$estimate, 2)
p_no_IL4   <- round(cor_no_IL4$p.value, 3)

cat("SENSITIVITY: Genomic PC1-LD50 excluding IL-4:\n")
cat("R =", r_no_IL4, "  p =", p_no_IL4, "\n")
