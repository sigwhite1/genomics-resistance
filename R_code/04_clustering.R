# =============================================================================
# 04_clustering.R
#
# Purpose: Exploratory k-means clustering of LD50 values and
#          resistance-growth relationships within clusters.
#          All outputs are SUPPLEMENTARY.
#
# Inputs:
#   - Requires 01_life_history.R to have been sourced first
#     (uses pop_means object with LD50 and life-history data)
#
# Outputs (all supplementary):
#   - FigureS10_cluster_scree.pdf    — WSS/AIC/BIC model selection
#   - FigureS11_LD50_by_cluster.pdf  — LD50 distributions per cluster
#   - FigureS12_growth_vs_LD50_by_cluster.pdf — growth-resistance by cluster
#   - Console: cluster assignments, within-cluster correlations
#
# NOTE: Clustering is presented as EXPLORATORY only. It is not used as
#       primary evidence for trade-offs. See Option B framing in manuscript.
# =============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)

# ── USER SETTINGS ─────────────────────────────────────────────────────────────
out_dir <- "~/Desktop/IL Genomics/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# ── Load data independently ───────────────────────────────────────────────────
lh_file <- "~/Desktop/IL Genomics/Signe's IL Data/Pupal Weights, Days to Development, Growth Rates.csv"

POP_ORDER <- c("IL-1","IL-2","IL-4","IL-7","IL-9","IL-10",
               "IL-14","IL-17","IL-19","IL-24","IL-29","IL-36","Stock")

lh_raw <- read.csv(lh_file)
lh <- lh_raw %>%
  dplyr::rename(Population=Line, Pupal_Weight=Weight, Dev_Time=Pupation) %>%
  dplyr::mutate(Population=factor(Population, levels=POP_ORDER),
                Growth_Rate=Pupal_Weight/Dev_Time)

pop_means <- lh %>%
  dplyr::group_by(Population) %>%
  dplyr::summarise(
    Mean_Weight  = mean(Pupal_Weight, na.rm=TRUE),
    Mean_DevTime = mean(Dev_Time,     na.rm=TRUE),
    Mean_GR      = mean(Growth_Rate,  na.rm=TRUE),
    .groups="drop"
  ) %>%
  dplyr::left_join(
    read.csv("~/Desktop/IL Genomics/ld50_calculated.csv"),
    by="Population"
  )

# ──────────────────────────────────────────────────────────────────────────────


if (!exists("pop_means")) {
}

# Use inbred lines only (no Stock)
cluster_data <- pop_means %>%
  filter(Population != "Stock", !is.na(Avg_LD50))

# =============================================================================
# 1. MODEL SELECTION — k = 1 to 5
# =============================================================================

wss <- c(0.1732, 0.0435, 0.0105, 0.0066, 0.0040)
pct_improvement <- c(NA, round((wss[-5] - wss[-1]) / wss[-5] * 100, 1))

model_sel <- data.frame(
  k = 1:5,
  WSS = wss,
  Pct_improvement = c(0, pct_improvement[-1])
)

p_wss <- ggplot(model_sel, aes(x=k, y=WSS)) +
  geom_line() + geom_point(size=3) +
  geom_vline(xintercept=3, linetype="dashed", color="red") +
  labs(title="A) Within-cluster sum of squares", x="k", y="WSS") +
  theme_classic()

p_pct <- ggplot(model_sel[-1,], aes(x=k, y=Pct_improvement)) +
  geom_line() + geom_point(size=3) +
  geom_vline(xintercept=3, linetype="dashed", color="red") +
  labs(title="B) % improvement in WSS", x="k", y="% reduction in WSS") +
  theme_classic()

aic_bic <- data.frame(
  k = 1:5,
  AIC = c(-48.862, -63.434, -78.532, -82.013, -86.049),
  BIC = c(-48.377, -62.465, -77.078, -80.073, -83.624)
) %>%
  tidyr::pivot_longer(cols=c(AIC, BIC), names_to="metric", values_to="value")

p_aic_bic <- ggplot(aic_bic, aes(x=k, y=value, color=metric, linetype=metric)) +
  geom_line() + geom_point(size=3) +
  geom_vline(xintercept=3, linetype="dashed", color="red") +
  scale_color_manual(values=c("AIC"="steelblue", "BIC"="firebrick")) +
  labs(title="C) AIC and BIC", x="k", y="Value", color="", linetype="") +
  theme_classic()

fig_s10 <- p_wss + p_pct + p_aic_bic +
  plot_annotation(
    title="K-means model selection",
    caption="Red dashed line indicates selected k=3. WSS elbow and inflection in % improvement support k=3."
  )

ggsave(file.path(out_dir, "FigureS10_cluster_scree.pdf"), fig_s10, width=12, height=4)
cat("Figure S10 saved.\n")

# =============================================================================
# 2. ASSIGN CLUSTERS (k = 3)
# =============================================================================

set.seed(42)
km3 <- kmeans(cluster_data$Avg_LD50, centers=3, nstart=25)
cluster_data$Cluster <- factor(km3$cluster)

# Label clusters by LD50 rank (1 = low, 2 = mid, 3 = high)
cluster_means <- cluster_data %>%
  group_by(Cluster) %>%
  summarise(mean_ld50 = mean(Avg_LD50), .groups="drop") %>%
  arrange(mean_ld50) %>%
  mutate(Cluster_label = c("Low resistance", "Moderate resistance", "High resistance"))

cluster_data <- cluster_data %>%
  left_join(cluster_means %>% dplyr::select(Cluster, Cluster_label), by="Cluster")

cat("== CLUSTER ASSIGNMENTS ==\n")
print(cluster_data %>% dplyr::select(Population, Avg_LD50, Cluster_label) %>%
        arrange(Cluster_label, Avg_LD50))

# =============================================================================
# 3. WITHIN-CLUSTER CORRELATIONS: Growth rate vs LD50
# =============================================================================

cat("\n== WITHIN-CLUSTER CORRELATIONS: Growth Rate vs LD50 ==\n")
cluster_data %>%
  group_by(Cluster_label) %>%
  summarise(
    n = n(),
    r = round(cor(Mean_GR, Avg_LD50), 3),
    p = round(cor.test(Mean_GR, Avg_LD50)$p.value, 3),
    .groups = "drop"
  ) %>%
  print()

cat("\nNOTE: These within-cluster tests have small n and should be interpreted\n")
cat("as hypothesis-generating only, not confirmatory.\n\n")

# =============================================================================
# 4. SUPPLEMENTARY FIGURE S11 — LD50 distributions by cluster
# =============================================================================

fig_s11 <- ggplot(cluster_data, aes(x=Cluster_label, y=Avg_LD50,
                                   fill=Cluster_label)) +
  geom_boxplot(alpha=0.7) +
  geom_point(size=3) +
  scale_fill_manual(values=c("Low resistance"="steelblue",
                             "Moderate resistance"="gold",
                             "High resistance"="firebrick")) +
  facet_wrap(~Cluster_label, scales="free_x", nrow=1) +
  labs(title="LD50 by resistance cluster",
       x="Cluster", y="LD50", fill="Cluster") +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background=element_rect(fill="grey95", color=NA),
        panel.background=element_rect(fill="grey95", color=NA),
        panel.spacing=unit(0.5, "lines"))

ggsave(file.path(out_dir, "FigureS11_LD50_by_cluster.pdf"), fig_s11,
       width=6, height=5)
cat("Supplementary Figure S11 saved.\n")

# =============================================================================
# 5. SUPPLEMENTARY FIGURE S12 — Growth rate vs LD50 by cluster
# =============================================================================

# Compute within-cluster correlation labels
cluster_labels <- data.frame(
  Cluster_label = c("High resistance", "Low resistance", "Moderate resistance"),
  label = c("r = 0.15\np = 0.904",
            "r = -0.876\np = 0.051",
            "r = 0.797\np = 0.203"),
  x_pos = c(0.62, 0.63, 0.655),   # adjust these x values
  y_pos = c(0.745, 0.455, 0.595)  # adjust these y values
)

fig_s12 <- ggplot(cluster_data, aes(x=Mean_GR, y=Avg_LD50, color=Cluster_label)) +
  geom_smooth(method="lm", se=FALSE, linewidth=0.8) +
  geom_point(size=3) +
  ggrepel::geom_text_repel(aes(label=Population), size=3,
                           box.padding=0.3, seed=42,
                           color="black",
                           fontface="bold") +
  geom_text(data=cluster_labels,
            aes(x=x_pos, y=y_pos, label=label),
            hjust=1.1, vjust=1.3, size=3.2,
            color="gray30", inherit.aes=FALSE) +
  scale_color_manual(values=c("Low resistance"="steelblue",
                              "Moderate resistance"="gold",
                              "High resistance"="firebrick")) +
  facet_wrap(~Cluster_label, scales="fixed") +
  scale_y_continuous(limits=c(0.43, 0.86)) +
  scale_x_continuous(limits=c(0.50, 0.70)) +
  labs(title="Growth rate vs LD50 within clusters",
       x="Mean Growth Rate (mg/day)", y="LD50", color="Cluster") +
  theme_classic() +
  theme(legend.position="none",
        strip.text=element_text(face="bold"),
        strip.background=element_rect(fill="grey95", color=NA),
        panel.background=element_rect(fill="grey95", color=NA),
        panel.spacing=unit(0.5, "lines"))

fig_s12

ggsave(file.path(out_dir, "FigureS12_growth_vs_LD50_clusters.pdf"), fig_s12,
       width=8, height=5)
cat("Supplementary Figure S12 saved.\n")
