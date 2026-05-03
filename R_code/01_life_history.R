# =============================================================================
# 01_life_history.R
#
# Purpose: Life-history PCA and figures
#
# Inputs:
#   - Pupal_Weights__Days_to_Development__Growth_Rates.csv
#   - Infection_Growth.csv
#
# Outputs:
#   - Figure_1_boxplots.pdf         — pupal weight, dev time, growth rate
#   - FigureS1_weight_vs_devtime.pdf — correlation between weight and dev time
#   - Figure_5_LH_PC1_vs_LD50.pdf  — life-history PC1 vs LD50 
#   - FigureS9_LD50_vs_traits.pdf  — LD50 vs individual life-history traits
#   - Console: PC1 loadings, ANOVA summaries, population means table
#
# NOTE: The headline PC1-LD50 result (genomic PC1, R = -0.63) is in
#       02_genomic_pca.R, not here. This script covers life-history only.
# =============================================================================

library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(emmeans)
library(multcomp)
library(MASS)
library(tibble)

# ── USER SETTINGS ─────────────────────────────────────────────────────────────
lh_file  <- "~/Desktop/Pupal Weights, Days to Development, Growth Rates.csv"
inf_file <- "~/Desktop/Infection_Growth.csv"
out_dir  <- "~/Desktop/figures"
# ──────────────────────────────────────────────────────────────────────────────

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Population order used consistently throughout all figures
POP_ORDER <- c("IL-1","IL-2","IL-4","IL-7","IL-9","IL-10",
               "IL-14","IL-17","IL-19","IL-24","IL-29","IL-36","Stock")

# Colour palette
POP_COLORS <- c(
  "Stock" = "white",   "IL-1"  = "#D73027", "IL-2"  = "#FC8D59",
  "IL-4"  = "#EFC000", "IL-7"  = "#91BFDB", "IL-9"  = "#4575B4",
  "IL-10" = "#313695", "IL-14" = "#8E44AD", "IL-17" = "#D73092",
  "IL-19" = "#666666", "IL-24" = "#1D91C0", "IL-29" = "#8C510A",
  "IL-36" = "#A6DBA0"
)

library(MASS)

bin_file <- "~/Desktop/Binomial Mortality Data.csv"
bin_raw  <- read.csv(bin_file)

populations <- c("IL-1","IL-2","IL-4","IL-7","IL-9","IL-10",
                 "IL-14","IL-17","IL-19","IL-24","IL-29","IL-36","Stock")

ld50_means <- data.frame(Population=populations, Avg_LD50=NA_real_, SE=NA_real_)

for (pop in populations) {
  col <- gsub("-", ".", pop)
  df  <- data.frame(infected=bin_raw[[col]], dose=bin_raw$Dose)
  df  <- df[df$dose > 0 & !is.na(df$infected), ]
  fit      <- glm(infected ~ dose, data=df, family=binomial)
  ld50_est <- dose.p(fit, p=0.5)
  ld50_means$Avg_LD50[ld50_means$Population == pop] <- as.numeric(ld50_est)
  ld50_means$SE[ld50_means$Population == pop]        <- as.numeric(attr(ld50_est, "SE"))
}

ld50_means$Population <- factor(ld50_means$Population, levels=POP_ORDER)
write.csv(ld50_means, "~/Desktop/IL Genomics/ld50_calculated.csv", row.names = FALSE)

# =============================================================================
# 1. LOAD AND CLEAN DATA
# =============================================================================

lh_raw <- read.csv(lh_file)

lh <- lh_raw %>%
  rename(Population = Line, Pupal_Weight = Weight, Dev_Time = Pupation) %>%
  mutate(
    Population  = factor(Population, levels = POP_ORDER),
    Growth_Rate = Pupal_Weight / Dev_Time
  )

cat("Life-history data loaded:", nrow(lh), "rows,",
    n_distinct(lh$Population), "populations\n")

# =============================================================================
# 2. FIGURE 1 — Boxplots of life-history traits by population
# =============================================================================

tukey_letters <- function(df, trait_col) {
  formula    <- as.formula(paste(trait_col, "~ Population"))
  model      <- aov(formula, data = df)
  cld_result <- multcomp::cld(emmeans(model, ~ Population),
                               Letters = letters, adjust = "tukey")
  cld_df <- as.data.frame(cld_result)
  if (!"Population" %in% colnames(cld_df) || !".group" %in% colnames(cld_df)) {
    colnames(cld_df)[1]         <- "Population"
    colnames(cld_df)[ncol(cld_df)] <- ".group"
  }
  cld_df %>% dplyr::select(Population, .group) %>% mutate(.group = trimws(.group))
}

pw_letters  <- tukey_letters(lh, "Pupal_Weight")
dev_letters <- tukey_letters(lh, "Dev_Time")
gr_letters  <- tukey_letters(lh, "Growth_Rate")

boxplot_theme <- theme_minimal() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.position = "none")

make_boxplot <- function(df, y_col, y_lab, title, letters_df) {
  y_max <- max(df[[y_col]], na.rm = TRUE)
  ggplot(df, aes(x = Population, y = .data[[y_col]], fill = Population)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.5, color = "grey30") +
    geom_text(data = letters_df,
              aes(x = Population, y = y_max * 1.05, label = .group),
              inherit.aes = FALSE, size = 3) +
    scale_fill_manual(values = POP_COLORS) +
    scale_x_discrete(limits = POP_ORDER) +
    labs(title = title, x = "Population", y = y_lab) +
    boxplot_theme
}

p1a <- make_boxplot(lh, "Pupal_Weight", "Pupal Weight (mg)",
                    "A) Pupal Weight by Population", pw_letters)
p1b <- make_boxplot(lh, "Dev_Time", "Days to Development",
                    "B) Development Time by Population", dev_letters)
p1c <- make_boxplot(lh, "Growth_Rate", "Growth Rate (mg/day)",
                    "C) Growth Rate by Population", gr_letters)

fig1 <- p1a / p1b / p1c
ggsave(file.path(out_dir, "Figure_1_boxplots.pdf"), fig1, width=8, height=14)
cat("Figure 1 saved.\n")

cat("\n-- ANOVA: Pupal Weight --\n")
print(summary(aov(Pupal_Weight ~ Population, data = lh)))
cat("\n-- ANOVA: Development Time --\n")
print(summary(aov(Dev_Time ~ Population, data = lh)))
cat("\n-- ANOVA: Growth Rate --\n")
print(summary(aov(Growth_Rate ~ Population, data = lh)))

# =============================================================================
# 3. POPULATION-LEVEL MEANS
# =============================================================================

pop_means <- lh %>%
  group_by(Population) %>%
  summarise(
    Mean_Weight  = mean(Pupal_Weight, na.rm = TRUE),
    SE_Weight    = sd(Pupal_Weight,   na.rm = TRUE) / sqrt(n()),
    Mean_DevTime = mean(Dev_Time,     na.rm = TRUE),
    SE_DevTime   = sd(Dev_Time,       na.rm = TRUE) / sqrt(n()),
    Mean_GR      = mean(Growth_Rate,  na.rm = TRUE),
    SE_GR        = sd(Growth_Rate,    na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  left_join(ld50_means, by = "Population")

# =============================================================================
# 4. FIGURE S1 — Pupal weight vs. development time correlation
# =============================================================================

cor_fig2 <- cor.test(pop_means$Mean_Weight, pop_means$Mean_DevTime)
r_fig2   <- round(cor_fig2$estimate, 2)
p_fig2   <- signif(cor_fig2$p.value, 3)

figs1 <- ggplot(pop_means, aes(x=Mean_DevTime, y=Mean_Weight)) +
  geom_smooth(method="lm", color="black", linewidth=0.8, se=TRUE, fill="grey85") +
  geom_errorbar(data=filter(pop_means, Population != "Stock"),
                aes(ymin=Mean_Weight-SE_Weight, ymax=Mean_Weight+SE_Weight,
                    color=Population), width=0) +
  geom_errorbarh(data=filter(pop_means, Population != "Stock"),
                 aes(xmin=Mean_DevTime-SE_DevTime, xmax=Mean_DevTime+SE_DevTime,
                     color=Population), height=0) +
  geom_errorbar(data=filter(pop_means, Population == "Stock"),
                aes(ymin=Mean_Weight-SE_Weight, ymax=Mean_Weight+SE_Weight),
                color="black", width=0) +
  geom_errorbarh(data=filter(pop_means, Population == "Stock"),
                 aes(xmin=Mean_DevTime-SE_DevTime, xmax=Mean_DevTime+SE_DevTime),
                 color="black", height=0) +
  geom_point(aes(fill=Population), size=3.5, shape=21, stroke=0.8, color="black") +
  scale_fill_manual(values=POP_COLORS, name="Population") +
  scale_color_manual(values=POP_COLORS, guide="none") +
  annotate("text", x=min(pop_means$Mean_DevTime),
           y=max(pop_means$Mean_Weight)*1.02,
           label=paste0("r = ", r_fig2, "\np = ", p_fig2), hjust=0, size=4) +
  labs(title="Correlation: Mean Pupal Weight vs. Mean Development Time",
       x="Mean Days to Development", y="Mean Pupal Weight (mg)",
       color="Population") +
  theme_minimal() +
  theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.5))

ggsave(file.path(out_dir, "FigureS1_weight_vs_devtime.pdf"), fig2, width=8, height=6)
cat("Figure S1 saved.  r =", r_fig2, "  p =", p_fig2, "\n")

# =============================================================================
# 5. LIFE-HISTORY PCA
# =============================================================================

pca_input <- pop_means %>%
  filter(Population != "Stock") %>%
  dplyr::select(Population, Mean_Weight, Mean_DevTime, Mean_GR) %>%
  column_to_rownames("Population")

lh_pca <- prcomp(pca_input, scale. = TRUE)

cat("\n== LIFE-HISTORY PCA ==\n")
pct_var <- round(summary(lh_pca)$importance[2, ] * 100, 1)
cat("Variance explained:\n"); print(pct_var)
cat("\nPC loadings:\n"); print(round(lh_pca$rotation, 3))
cat("\nAll loadings on PC1 are negative — PC1 reflects the FAST life-history axis.\n")
cat("Higher PC1 = faster developing, lighter, higher growth rate.\n\n")

pc_scores <- as.data.frame(lh_pca$x) %>%
  tibble::rownames_to_column("Population")

pop_means <- pop_means %>%
  left_join(pc_scores %>% dplyr::select(Population, PC1, PC2), by="Population")

# =============================================================================
# 6. FIGURE 5 — Life-history PC1 vs. LD50
# =============================================================================

fig_5 <- ggplot(fig8_data, aes(x=PC1, y=Avg_LD50)) +
  geom_smooth(method="lm", color="black", linetype="solid",
              se=TRUE, fill="grey85", linewidth=0.9) +
  geom_point(aes(fill=Avg_LD50), size=3.5, shape=21, stroke=0.5) +
  scale_fill_gradientn(colors=c("steelblue","gold","firebrick"),
                       guide="none") +
  geom_text_repel(aes(label=Population), size=3, box.padding=0.4,
                  seed=42, fontface="italic") +
  annotate("text", x=Inf, y=-Inf,
           label=paste0("r = ", r_lh, ",  p = ", p_lh),
           hjust=1.1, vjust=-0.5, size=4) +
  labs(title="Life-history PC1 vs. LD50",
       x=paste0("Life-history PC1 (", pct_var["PC1"], "% variance)"),
       y="LD50") +
  theme_classic()

ggsave(file.path(out_dir, "Figure_5_LH_PC1_vs_LD50.pdf"), fig_s8, width=7, height=6)
cat("Figure 5 saved.\n")
# =============================================================================
# 7. SUPPLEMENTARY FIGURE S9 — LD50 vs. individual life-history traits
# =============================================================================

make_cor_plot <- function(df, x_col, x_lab, title) {
  ct      <- cor.test(df[[x_col]], df$Avg_LD50)
  r       <- round(ct$estimate, 3)
  p       <- signif(ct$p.value, 3)
  x_range <- range(df[[x_col]], na.rm=TRUE)
  y_range <- range(df$Avg_LD50,  na.rm=TRUE)
  ggplot(df, aes(x=.data[[x_col]], y=Avg_LD50)) +
    geom_smooth(method="lm", color="black", linetype="dashed",
                se=FALSE, linewidth=0.8) +
    geom_point(size=2.5, color="grey30") +
    annotate("text", x=x_range[2], y=y_range[1],
             label=paste0("r = ", r, "\np = ", p),
             hjust=1.25, vjust=0, size=3.5) +
    labs(title=title, x=x_lab, y="LD50") +
    theme_classic()
}

s9_data <- pop_means %>% filter(!is.na(Avg_LD50))

s9a <- make_cor_plot(s9_data, "Mean_GR",     "Growth Rate (mg/day)",       "A) LD50 vs. Growth Rate")
s9b <- make_cor_plot(s9_data, "Mean_Weight",  "Mean Pupal Weight (mg)",     "B) LD50 vs. Pupal Weight")
s9c <- make_cor_plot(s9_data, "Mean_DevTime", "Mean Development Time (days)","C) LD50 vs. Development Time")

fig_s9 <- s9a / s9b / s9c +
  plot_annotation(title="LD50 vs. Life-History Traits")

fig_s9

ggsave(file.path(out_dir, "FigureS9_LD50_vs_traits.pdf"), fig_s9, width=6, height=12)
cat("Supplementary Figure S9 saved.\n")

# =============================================================================
# 8. SUMMARY TABLE
# =============================================================================

cat("\n== POPULATION MEANS SUMMARY ==\n")
print(
  pop_means %>%
    dplyr::select(Population, Mean_Weight, SE_Weight, Mean_DevTime, SE_DevTime,
                  Mean_GR, SE_GR, Avg_LD50, PC1) %>%
    mutate(across(where(is.numeric), ~round(.x, 4))) %>%
    arrange(Population),
  n = Inf
)

cat("\nAll figures saved to:", normalizePath(out_dir), "\n")
cat("Script complete.\n")

# =============================================================================
# INFECTION ASSAY GLMs (Supplementary Table S3)
# =============================================================================

inf_raw <- read.csv("~/Desktop/IL Genomics/Signe's IL Data/Infection_Growth.csv")
inf_raw$Log_Dose   <- log10(inf_raw$Dose)
inf_raw$Population <- factor(inf_raw$Line, levels=POP_ORDER)

# Binomial logistic GLM
glm_logit <- glm(
  cbind(Infection_Status, Total_Exposed - Infection_Status) ~
    Log_Dose * Population,
  family = binomial,
  data   = inf_raw
)

cat("\n── Supplementary Table S3: ANOVA — Binomial logistic GLM ──\n")
print(anova(glm_logit, test = "Chisq"))

# =============================================================================
# 9. SUPPLEMENTARY FIGURES S4 - S6 — Dose-response visualisations
# =============================================================================

library(tidyr)

# Predicted infection probabilities across dose range (fine grid)
dose_seq  <- 10^seq(log10(min(inf_raw$Dose)), log10(max(inf_raw$Dose)), length.out=200)
pred_grid <- expand.grid(
  Dose       = dose_seq,
  Population = factor(levels(inf_raw$Population), levels=POP_ORDER)
)
pred_grid$Log_Dose <- log10(pred_grid$Dose)

pred_out      <- predict(glm_logit, newdata=pred_grid, type="link", se.fit=TRUE)
pred_grid$fit <- pred_out$fit
pred_grid$se  <- pred_out$se.fit
pred_grid$pred_prob  <- plogis(pred_grid$fit)
pred_grid$lower      <- plogis(pred_grid$fit - 1.96 * pred_grid$se)
pred_grid$upper      <- plogis(pred_grid$fit + 1.96 * pred_grid$se)
pred_grid$log_odds   <- pred_grid$fit   # log-odds = linear predictor

# Observed proportions
obs <- inf_raw %>%
  mutate(Observed = Infection_Status / Total_Exposed)

# ── Figure S4: All populations overlaid, proportion infected vs dose ──────────

fig_s4 <- ggplot() +
  geom_ribbon(data = pred_grid,
              aes(x = log10(Dose), ymin = lower, ymax = upper, fill = Population),
              alpha = 0.15) +
  geom_line(data = pred_grid,
            aes(x = log10(Dose), y = pred_prob, color = Population),
            linewidth = 0.9) +
  geom_point(data = obs,
             aes(x = log10(Dose), y = Observed, color = Population),
             size = 2.5) +
  scale_color_manual(values = POP_COLORS, breaks = POP_ORDER) +
  scale_fill_manual(values = POP_COLORS, breaks = POP_ORDER, guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title  = "Dose-Response Curve (Logistic GLM)",
    x      = expression(Log[10](Dose)),
    y      = "Proportion Infected",
    color  = "Population"
  ) +
  theme_classic() +
  theme(legend.position = "right")

fig_s4

ggsave(file.path(out_dir, "FigureS4_dose_response_overlaid.pdf"), fig_s4,
       width = 10, height = 6)


# ── Figure S5: Observed vs predicted, faceted by population ──────────────────

fig_s5 <- ggplot() +
  geom_ribbon(data = pred_grid,
              aes(x = log10(Dose), ymin = lower, ymax = upper, fill = Population),
              alpha = 0.2) +
  geom_line(data = pred_grid,
            aes(x = log10(Dose), y = pred_prob, color = Population),
            linewidth = 0.9) +
  geom_line(data = obs,
            aes(x = log10(Dose), y = Observed, group = Population),
            linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_point(data = obs,
             aes(x = log10(Dose), y = Observed),
             color = "red", size = 2) +
  # Manual legend for line types
  geom_line(data = pred_grid[1,],
            aes(x = log10(Dose), y = pred_prob, linetype = "Predicted"),
            color = "black") +
  geom_point(data = obs[1,],
             aes(x = log10(Dose), y = Observed, shape = "Observed"),
             color = "red") +
  scale_linetype_manual(name = "Type",
                        values = c("Predicted" = "solid")) +
  scale_shape_manual(name = "Type",
                     values = c("Observed" = 16)) +
  scale_fill_manual(values = POP_COLORS, guide = "none") +
  scale_color_manual(values = POP_COLORS, guide = "none") +
  facet_wrap(~ Population, ncol = 4) +
  labs(
    title = "Observed vs. Predicted Infection Rates by Population",
    x     = expression(Log[10](Dose)),
    y     = "Infection Proportion"
  ) +
  theme_classic() +
  theme(
    strip.text      = element_text(face = "italic"),
    legend.title    = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  guides(
    linetype = guide_legend(override.aes = list(color = "black", linewidth = 0.9)),
    shape    = guide_legend(override.aes = list(color = "red", size = 2))
  )
fig_s5

ggsave(file.path(out_dir, "FigureS5_obs_vs_pred.pdf"), fig_s5,
       width=9, height=6)
cat("Supplementary Figure S5 saved.\n")

# ── Figure S6: Log-odds vs log10(dose), all populations overlaid ─────────────

fig_s6 <- ggplot(pred_grid, aes(x=log10(Dose), y=log_odds, color=Population)) +
  geom_line(linewidth=0.9) +
  scale_color_manual(values=POP_COLORS) +
  labs(
    title = "Log-Odds vs Log10(Dose) by Population",
    x     = expression(Log[10](Dose)),
    y     = "Log-Odds of Infection",
    color = "Population"
  ) +
  theme_classic() +
  theme(legend.position="right")

fig_s6

ggsave(file.path(out_dir, "FigureS6_logodds_dose.pdf"), fig_s6,
       width=10, height=6)
cat("Supplementary Figure S6 saved.\n")
