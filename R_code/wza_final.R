# =============================================================================
# wza_final.R
#
# Final WZA + PCA-LD50 joint analysis under the pre-specified inclusion
# criteria and corrected normalization.
#
# INCLUSION CRITERIA (decided from window properties, not from results):
#   1. Full-size windows only (exactly 250 kb). Truncated chromosome-end
#      fragments are not comparable units.
#   2. At least 100 variants per window. Seven full-size windows fall below
#      this, clustering into contiguous variant-depleted regions on
#      NC_071299.1, NC_071306.1 and NC_071296.1 at 20-50x below the
#      genome-wide median (2077). WZA is unstable when computed from so few
#      markers, and such regions carry little variation to associate with
#      phenotype.
#   Retains 2164 of 2232 windows (97%).
#
# NORMALIZATION FIX:
#   Dispersion model fitted as a Gamma GLM with log link so predicted SDs are
#   guaranteed positive (the original lm() produced negative SDs, and
#   division by near-zero SDs produced |Z| up to 161). Mean absolute
#   deviation converted to SD via sigma = MAD * sqrt(pi/2).
# =============================================================================

.libPaths("~/R/library")
library(dplyr)
setwd("/global/scratch/users/pswhite/analysis")

MIN_SNPS    <- 100
WINDOW_SIZE <- 250000

cat("Loading WZA cache...\n")
f <- list.files("wza_cache", pattern = "*.txt$", full.names = TRUE)
d <- bind_rows(lapply(f, read.delim))

d$ep <- rank(d$p, ties.method = "average") / nrow(d)
d$z  <- qnorm(1 - d$ep)
d$He <- d$maf * (1 - d$maf)

w <- d %>%
  group_by(window) %>%
  summarise(n_SNPs = n(), ZW = sum(He * z) / sqrt(sum(He^2)), .groups = "drop") %>%
  filter(is.finite(ZW))

w$start <- as.numeric(sub(".*_([0-9]+)_[0-9]+\\.vcf$", "\\1", w$window))
w$end   <- as.numeric(sub(".*_[0-9]+_([0-9]+)\\.vcf$", "\\1", w$window))
w$width <- w$end - w$start

n_before <- nrow(w)
w <- w %>% filter(width == WINDOW_SIZE, n_SNPs >= MIN_SNPS)
cat("Windows before filtering:", n_before, " | after:", nrow(w),
    "(", round(100*nrow(w)/n_before, 1), "% )\n")

# --- Corrected normalization ---
mean_model <- lm(ZW ~ poly(n_SNPs, 2), data = w)
w$pred_mean <- predict(mean_model)
w$abs_resid <- abs(w$ZW - w$pred_mean) # per-window absolute residual (mean AD, not median AD)

sd_model <- glm(pmax(abs_resid, 1e-8) ~ poly(n_SNPs, 2),
                family = Gamma(link = "log"), data = w)
w$pred_SD <- predict(sd_model, type = "response") * sqrt(pi / 2)

w$ZW_norm  <- (w$ZW - w$pred_mean) / w$pred_SD
w$abs_ZW   <- abs(w$ZW_norm)
w$wza_rank <- rank(-w$abs_ZW, ties.method = "min")
w$clean_window <- gsub("^window_", "", w$window)

cat("Predicted SD range:", signif(range(w$pred_SD), 4),
    " | any <= 0:", any(w$pred_SD <= 0), "\n")
cat("max |Z_W norm|:", round(max(w$abs_ZW), 2), "\n")

# --- PCA-LD50 side, restricted to the same window set ---
snprelate <- read.csv("snprelate_1mb.csv")
ld50 <- read.csv("ld50_calculated.csv")
ld50 <- ld50[ld50$Population != "Stock", c("Population","Avg_LD50")]
colnames(ld50) <- c("sample_id","ld50")

m <- merge(snprelate, ld50, by = "sample_id")
rownames(m) <- m$sample_id; m$sample_id <- NULL
m$ld50 <- as.numeric(m$ld50)

wc <- setdiff(colnames(m), "ld50")
cv <- sapply(wc, function(col) cor.test(m$ld50, m[[col]])$estimate)
pv <- sapply(wc, function(col) cor.test(m$ld50, m[[col]])$p.value)

pve <- data.frame(clean_window = gsub("slidingwindow_window_", "", wc),
                  correlation = as.numeric(cv),
                  p_value = as.numeric(pv),
                  stringsAsFactors = FALSE)
pve$PVE <- pve$correlation^2
pve <- pve[!is.na(pve$PVE), ]

both <- merge(pve, w[, c("clean_window","ZW_norm","abs_ZW","wza_rank","n_SNPs")],
              by = "clean_window")
both$pve_rank <- rank(-both$PVE, ties.method = "min")

n <- nrow(both)
t1 <- ceiling(0.01*n); t2 <- ceiling(0.02*n)
cat("\nWindows in joint analysis:", n, " | top1% <=", t1, " | top2% <=", t2, "\n")
cat("Min PCA-LD50 p-value:", signif(min(both$p_value), 3),
    " | Bonferroni threshold:", signif(0.05/n, 3), "\n")

both$tier <- ifelse(both$pve_rank <= t1 & both$wza_rank <= t1, "Top 1% (both)",
             ifelse(both$pve_rank <= t2 & both$wza_rank <= t2, "Top 2% (both)", NA))

conv <- both %>% filter(!is.na(tier)) %>% arrange(pve_rank)

cat("\n=== FINAL CONVERGENT WINDOWS ===\n")
if (nrow(conv) == 0) {
  cat("No windows meet the joint criterion.\n")
} else {
  print(as.data.frame(conv[, c("clean_window","tier","correlation","PVE",
                               "pve_rank","ZW_norm","wza_rank","n_SNPs")]),
        row.names = FALSE)
}

write.csv(both, "joint_results_FINAL.csv", row.names = FALSE)
write.csv(conv, "convergent_windows_FINAL.csv", row.names = FALSE)
cat("\nSaved joint_results_FINAL.csv and convergent_windows_FINAL.csv\n")
