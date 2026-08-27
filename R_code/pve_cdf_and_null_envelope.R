# =============================================================================
# pve_cdf_and_null_envelope.R
#
# Purpose: (1) Rebuild the PVE cumulative distribution correctly, from the
#          verified 2,164-window joint_results_FINAL.csv (NOT from pve_cdf.csv,
#          which appears to have been built pre-inclusion-criteria from the
#          raw ~2,230-window `pve` object in wza_final.R and is therefore stale
#          -- see chat discussion).
#          (2) Answer Reviewer 2's implicit concern that a diffuse PVE curve
#          could also arise from pure noise (no true signal at all, n=12).
#          Do this by permuting LD50 across the 12 lines 1,000 times,
#          recomputing the PVE distribution each time using the identical
#          data-loading logic wza_final.R uses for the real analysis, and
#          building a null envelope. If the observed curve falls outside the
#          null envelope, diffuse-but-real signal is distinguished from
#          diffuse noise directly, rather than argued for.
#
# Prerequisites (Savio):
#   - joint_results_FINAL.csv   (from wza_final.R; source of the real PVE curve
#                                 and the authoritative 2,164-window set)
#   - snprelate_1mb.csv, ld50_calculated.csv  (same inputs wza_final.R uses to
#                                 build `pve`; reused here for permutation)
#
# Output:
#   - pve_cdf_with_null_envelope.csv  (window_rank, observed_cumulative_share,
#     null_mean, null_lower_2.5, null_upper_97.5)
# =============================================================================

.libPaths("~/R/library")
setwd("/global/scratch/users/pswhite/analysis")
set.seed(42)

# -----------------------------------------------------------------------
# 1. REAL (OBSERVED) PVE CDF -- from the verified, criteria-filtered set
# -----------------------------------------------------------------------

joint <- read.csv("joint_results_FINAL.csv", stringsAsFactors = FALSE)
cat("joint_results_FINAL.csv rows:", nrow(joint), "(expect 2164)\n")

obs_pve_sorted <- sort(joint$PVE, decreasing = TRUE)
obs_cum_share  <- cumsum(obs_pve_sorted) / sum(obs_pve_sorted)

n_win <- length(obs_cum_share)
cat("Observed: top 1% (", ceiling(0.01 * n_win), "windows) ->",
    round(100 * obs_cum_share[ceiling(0.01 * n_win)], 1), "% of total PVE\n")
cat("Observed: top 25% (", ceiling(0.25 * n_win), "windows) ->",
    round(100 * obs_cum_share[ceiling(0.25 * n_win)], 1), "% of total PVE\n")

# -----------------------------------------------------------------------
# 2. NULL ENVELOPE -- same data-loading logic as wza_final.R's `pve` block,
#    restricted to the same 2,164-window set, with LD50 permuted across lines
# -----------------------------------------------------------------------

snprelate <- read.csv("snprelate_1mb.csv", stringsAsFactors = FALSE)
ld50      <- read.csv("ld50_calculated.csv", stringsAsFactors = FALSE)
ld50 <- ld50[ld50$Population != "Stock", c("Population", "Avg_LD50")]
colnames(ld50) <- c("sample_id", "ld50")

m <- merge(snprelate, ld50, by = "sample_id")
rownames(m) <- m$sample_id
m$sample_id <- NULL
m$ld50 <- as.numeric(m$ld50)
cat("\nMerged sample rows:", nrow(m), "(expect 12)\n")

wc <- setdiff(colnames(m), "ld50")
clean_names_all <- gsub("slidingwindow_window_", "", wc)

# Restrict to exactly the windows retained in joint_results_FINAL.csv
keep <- clean_names_all %in% joint$clean_window
wc_final <- wc[keep]
cat("Window columns matched to the final 2,164-window set:", length(wc_final),
    "(expect 2164)\n")
if (length(wc_final) != nrow(joint)) {
  cat("WARNING: window-name matching did not fully align -- check clean_window",
      "formatting between snprelate_1mb.csv and joint_results_FINAL.csv before",
      "trusting the null envelope below.\n")
}

window_matrix <- as.matrix(m[, wc_final, drop = FALSE])
storage.mode(window_matrix) <- "numeric"

n_perm <- 1000
null_cum_matrix <- matrix(NA_real_, nrow = n_perm, ncol = length(wc_final))

for (i in seq_len(n_perm)) {
  ld50_perm <- sample(m$ld50)                       # shuffle LD50 across the 12 lines
  cv_perm   <- as.numeric(cor(window_matrix, ld50_perm))
  pve_perm  <- cv_perm^2
  pve_perm_sorted <- sort(pve_perm, decreasing = TRUE)
  null_cum_matrix[i, ] <- cumsum(pve_perm_sorted) / sum(pve_perm_sorted)
}

null_mean  <- colMeans(null_cum_matrix)
null_lower <- apply(null_cum_matrix, 2, quantile, probs = 0.025)
null_upper <- apply(null_cum_matrix, 2, quantile, probs = 0.975)

# -----------------------------------------------------------------------
# 3. WHERE DOES THE OBSERVED CURVE EXIT THE NULL ENVELOPE, AND WHICH WAY?
#    Do NOT presume the direction. With n=12, pure noise across 2,164
#    windows can itself produce a few by-chance high-r^2 windows (a
#    multiple-testing/winner's-curse effect), which could make the NULL
#    curve look concentrated at the very top ranks. So:
#      observed ABOVE null envelope  -> MORE concentrated than noise
#                                        (argues against polygenicity)
#      observed BELOW null envelope  -> LESS concentrated / more evenly
#                                        spread than noise (consistent
#                                        with many small, fairly uniform
#                                        true effects)
#      observed WITHIN null envelope -> not distinguishable from pure
#                                        noise at this sample size
# -----------------------------------------------------------------------

above_upper <- which(obs_cum_share > null_upper)
below_lower <- which(obs_cum_share < null_lower)

first_above <- if (length(above_upper) > 0) min(above_upper) else NA
first_below <- if (length(below_lower) > 0) min(below_lower) else NA

cat("\nObserved curve exceeds null upper (97.5th pctile) starting at rank:",
    first_above,
    if (!is.na(first_above)) paste0(" (", round(100*first_above/n_win,1), "% of windows)") else "(never)",
    "\n")
cat("Observed curve falls below null lower (2.5th pctile) starting at rank:",
    first_below,
    if (!is.na(first_below)) paste0(" (", round(100*first_below/n_win,1), "% of windows)") else "(never)",
    "\n")
cat("Fraction of ranks where observed sits strictly inside the null envelope:",
    round(mean(obs_cum_share >= null_lower & obs_cum_share <= null_upper), 3), "\n")
cat("\nInterpret using the comment block above this section -- direction matters",
    "and should not be assumed before reading these numbers.\n")

# -----------------------------------------------------------------------
# 4. WRITE OUTPUT FOR PLOTTING
# -----------------------------------------------------------------------

out <- data.frame(
  window_rank              = seq_len(n_win),
  window_fraction          = seq_len(n_win) / n_win,
  observed_cumulative_PVE  = obs_cum_share,
  null_mean_cumulative_PVE = null_mean,
  null_lower_2.5           = null_lower,
  null_upper_97.5          = null_upper
)

write.csv(out, "pve_cdf_with_null_envelope.csv", row.names = FALSE)
cat("\nWrote pve_cdf_with_null_envelope.csv (", nrow(out), "rows )\n")
