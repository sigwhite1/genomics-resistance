# =============================================================================
# 03_manhattan_wza.R
#
# Purpose: Generate Manhattan plots, WZA scores, QQ plots, and top windows
#          table. Runs on Savio (requires pre-computed data).
#
# Prerequisites:
#   - wza_cache/ folder populated (run parse_vcf_final.R on Savio first)
#   - snprelate_1mb.csv generated (run parse_vcf_final.R on Savio first)
#
# Outputs:
#   - /global/scratch/users/analysis/figures.pdf
#     Page 1: Combined Manhattan plot (PCA-LD50 + WZA)
#     Page 2: QQ plots (PCA-LD50 + WZA)
#     Page 3: Top 10 genomic windows table
# =============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)
library(gridExtra)
library(grid)


manhattan_plot_data <- readRDS("~/Desktop/IL Genomics/manhattan_plot_data.rds")
pca_plot_data       <- readRDS("~/Desktop/IL Genomics/pca_plot_data.rds")
qq_wza              <- readRDS("~/Desktop/IL Genomics/qq_wza.rds")
qq_pca_data         <- readRDS("~/Desktop/IL Genomics/qq_pca_data.rds")
top10               <- readRDS("~/Desktop/IL Genomics/top10.rds")

# ── USER SETTINGS ─────────────────────────────────────────────────────────────
analysis_dir <- "/global/scratch/users/analysis"
out_pdf      <- file.path(analysis_dir, "figures.pdf")
# ──────────────────────────────────────────────────────────────────────────────

# LD50 values calculated in 01_life_history.R via MASS::dose.p
# (do NOT hardcode here — run 01_life_history.R first, then scp the CSV to Savio)
ld50_file <- "/global/scratch/users/analysis/ld50_calculated.csv"
if (!file.exists(ld50_file)) stop("ld50_calculated.csv not found. Upload from local Mac first.")
ld50_table  <- read.csv(ld50_file)
ld50_table  <- ld50_table[ld50_table$Population != "Stock", ]   # inbred lines only
sample_ids  <- ld50_table$Population
ld50_values <- ld50_table$Avg_LD50

# =============================================================================
# 1. LOAD AND PROCESS WZA CACHE DATA
# =============================================================================

cat("Loading WZA cache...\n")
cached_files <- list.files(file.path(analysis_dir, "wza_cache"),
                           pattern = "*.txt$", full.names = TRUE)
snp_data <- bind_rows(lapply(cached_files, read.delim))
cat("Loaded", nrow(snp_data), "SNPs from", length(cached_files), "windows\n")

# Compute empirical p-values and WZA components
snp_data$empirical_p <- rank(snp_data$p, ties.method="average") / nrow(snp_data)
snp_data$z           <- qnorm(1 - snp_data$empirical_p)
snp_data$He          <- snp_data$maf * (1 - snp_data$maf)

# WZA per window
wza_results <- snp_data %>%
  group_by(window) %>%
  summarise(
    n_SNPs    = n(),
    He_sq_sum = sum(He^2),
    ZW        = sum(He * z) / sqrt(He_sq_sum)
  ) %>%
  filter(is.finite(ZW))

# Normalize ZW by SNP count
mean_model <- lm(ZW ~ poly(n_SNPs, 2), data=wza_results)
sd_model   <- lm(abs(ZW - predict(mean_model)) ~ poly(n_SNPs, 2), data=wza_results)

wza_results$meanZW  <- predict(mean_model)
wza_results$sdZW    <- predict(sd_model)
wza_results$ZW_norm <- (wza_results$ZW - wza_results$meanZW) / wza_results$sdZW
wza_results$p_value <- 2 * (1 - pnorm(abs(wza_results$ZW_norm)))
wza_results$log10p  <- -log10(wza_results$p_value)

cat("WZA computed for", nrow(wza_results), "windows\n")

# =============================================================================
# 2. LOAD AND PROCESS SNPRELATE PCA DATA
# =============================================================================

cat("Loading SNPRelate PCA results...\n")
ordered <- read.csv(file.path(analysis_dir, "snprelate_1mb.csv"))
ordered$sample_id <- sample_ids

ld50_table <- data.frame(sample_id=sample_ids, ld50=ld50_values)
merged_data <- merge(ordered, ld50_table, by="sample_id")
merged_data <- merged_data %>% dplyr::select(sample_id, ld50, everything())
rownames(merged_data) <- merged_data$sample_id
merged_data <- merged_data[, -1]
merged_data$ld50 <- as.numeric(merged_data$ld50)

# Correlate each window PC1 with LD50
cor_values <- numeric(ncol(merged_data) - 1)
p_values   <- numeric(ncol(merged_data) - 1)

for (i in 1:(ncol(merged_data) - 1)) {
  test_result  <- cor.test(merged_data$ld50, merged_data[, i+1], method="pearson")
  cor_values[i] <- test_result$estimate
  p_values[i]   <- test_result$p.value
}

cor_results <- data.frame(
  window      = colnames(merged_data)[-1],
  correlation = cor_values,
  p_value     = p_values
)
cor_results <- cor_results[-ncol(merged_data), ]
cor_results$fdr_adj <- p.adjust(cor_results$p_value, method="fdr")
cor_results <- na.omit(cor_results)

fdr_threshold_logP <- -log10(0.05 / nrow(cor_results))
cat("Best window R:", round(max(abs(cor_results$correlation)), 3),
    "  p:", round(min(cor_results$p_value), 4), "\n")

# =============================================================================
# 3. BUILD MANHATTAN PLOT DATA
# =============================================================================

# WZA Manhattan
manhattan_data <- wza_results %>%
  mutate(
    chrom       = sub("window_((NC|NW|CM)_[0-9]+\\.[0-9]+)_.*", "\\1", window),
    position    = as.numeric(sub(".*_([0-9]+)_([0-9]+)\\.vcf$", "\\1", window)),
    end_position= as.numeric(sub(".*_([0-9]+)_([0-9]+)\\.vcf$", "\\2", window)),
    chrom_index = as.numeric(factor(chrom, levels=unique(chrom))),
    color       = as.factor(chrom_index %% 2),
    logP        = -log10(p_value)
  )

fdr_threshold_wza <- -log10(0.05 / nrow(manhattan_data))

max_window_pos <- manhattan_data %>%
  group_by(chrom) %>%
  summarise(max_pos=max(end_position, na.rm=TRUE)) %>%
  ungroup()
global_max_pos <- max(max_window_pos$max_pos)

padding_wza <- manhattan_data %>%
  group_by(chrom) %>%
  summarise() %>%
  mutate(position=global_max_pos, logP=NA, color="gray",
         window=NA, end_position=global_max_pos,
         chrom_index=as.numeric(factor(chrom, levels=unique(manhattan_data$chrom))))

manhattan_plot_data <- bind_rows(manhattan_data, padding_wza) %>%
  mutate(color=as.factor(ifelse(is.na(color), "padding", as.character(color))))

p_wza <- ggplot(manhattan_plot_data, aes(x=position, y=logP, color=color)) +
  geom_point(alpha=0.75, size=1.5, na.rm=TRUE) +
  geom_hline(yintercept=fdr_threshold_wza, linetype="dashed", color="red", size=1) +
  facet_grid(~chrom, scales="fixed", space="fixed") +
  scale_color_manual(values=c("0"="gray30","1"="gray60","padding"=NA), guide="none") +
  labs(title="WZA scores by chromosome",
       x="Chromosome", y="-log10(WZA P-value)") +
  theme_minimal() +
  theme(strip.text=element_text(size=6, angle=45, hjust=1, vjust=1),
        axis.text.x=element_blank(), axis.ticks.x=element_blank())

# PCA Manhattan
pca_manhattan_data <- cor_results %>%
  mutate(
    chrom       = sub("slidingwindow_window_((NC|NW)_[0-9]+\\.[0-9]+)_.*", "\\1", window),
    position    = as.numeric(sub(".*_([0-9]+)_([0-9]+)\\.vcf$", "\\1", window)),
    end_position= as.numeric(sub(".*_([0-9]+)_([0-9]+)\\.vcf$", "\\2", window)),
    chrom_index = as.numeric(factor(chrom, levels=unique(chrom))),
    color       = as.factor(chrom_index %% 2),
    logP        = -log10(p_value)
  )

max_window_pos2 <- pca_manhattan_data %>%
  group_by(chrom) %>%
  summarise(max_pos=max(end_position, na.rm=TRUE), .groups="drop")
global_max_pos2 <- max(max_window_pos2$max_pos)

padding_pca <- pca_manhattan_data %>%
  distinct(chrom, chrom_index) %>%
  mutate(position=global_max_pos2, end_position=global_max_pos2,
         logP=NA, color="padding", window=NA)

pca_plot_data <- bind_rows(pca_manhattan_data, padding_pca) %>%
  mutate(color=factor(color, levels=c("0","1","padding")))

p_pca <- ggplot(pca_plot_data, aes(x=position, y=logP, color=color)) +
  geom_point(alpha=0.75, size=1.5, na.rm=TRUE) +
  geom_hline(yintercept=fdr_threshold_logP, linetype="dashed", color="red", size=1) +
  facet_grid(~chrom, scales="fixed", space="fixed") +
  scale_color_manual(values=c("0"="gray30","1"="gray60","padding"=NA), guide="none") +
  labs(title="Manhattan Plot of PCA-LD50 Correlation",
       x="Chromosome", y=expression(-log[10](italic(P)))) +
  theme_minimal() +
  theme(strip.text=element_text(size=6, angle=45, hjust=1, vjust=1),
        axis.text.x=element_blank(), axis.ticks.x=element_blank())

# =============================================================================
# 4. QQ PLOTS
# =============================================================================

qq_wza <- wza_results %>%
  filter(is.finite(p_value), !is.na(p_value)) %>%
  arrange(p_value) %>%
  mutate(observed=-log10(p_value), expected=-log10(ppoints(n())))

p_qq_wza <- ggplot(qq_wza, aes(x=expected, y=observed)) +
  geom_point(size=1.5, alpha=0.7, color="blue") +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red") +
  theme_minimal() +
  labs(title="WZA QQ Plot",
       x="Expected -log10(P)", y="Observed -log10(P)")

observed_pvals_log <- -log10(sort(cor_results$p_value))
expected_pvals     <- -log10(ppoints(length(observed_pvals_log)))
qq_pca_data        <- data.frame(Expected=expected_pvals, Observed=observed_pvals_log)

p_qq_pca <- ggplot(qq_pca_data, aes(x=Expected, y=Observed)) +
  geom_point(alpha=0.7, size=2, color="blue") +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
  theme_minimal() +
  labs(title="PCA-LD50 QQ Plot",
       x="Expected -log10(P)", y="Observed -log10(P)")

# =============================================================================
# 5. TOP 10 WINDOWS TABLE
# =============================================================================

wza_results$wza_fdr  <- p.adjust(wza_results$p_value, method="fdr")
wza_results$wza_rank <- rank(wza_results$p_value, ties.method="min")
cor_results$pca_fdr  <- p.adjust(cor_results$p_value, method="fdr")
cor_results$pca_rank <- rank(cor_results$p_value, ties.method="min")

wza_results$clean_window <- gsub("^window_", "", wza_results$window)
cor_results$clean_window <- gsub("slidingwindow_window_", "", cor_results$window)

combined_df <- merge(
  wza_results[, c("clean_window","wza_fdr","wza_rank")],
  cor_results[, c("clean_window","pca_fdr","pca_rank")],
  by="clean_window"
)

top10 <- combined_df[order(combined_df$wza_rank), ][1:10, ]
colnames(top10) <- c("Window","WZA FDR","WZA Rank","PCA FDR","PCA Rank")
top10$`WZA FDR` <- signif(top10$`WZA FDR`, 3)
top10$`PCA FDR` <- signif(top10$`PCA FDR`, 3)

cat("\nTop 10 genomic windows by WZA:\n")
print(top10)

# =============================================================================
# 6. SAVE ALL TO PDF
# =============================================================================

cat("\nSaving figures to", out_pdf, "...\n")
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
cat("Done! Figures saved to", out_pdf, "\n")
