# =============================================================================
# 05_CDS_vs_genomewide_comparison.R
#
# Complete CDS-vs-genome-wide robustness check (Reviewer 3, CDS-filtering
# comment). Two parts, run in two different environments:
#
#   PART 1 (Savio):  genome-wide glPca -- produces genomewide_pca_scores.csv
#                     (needed for Figure S13, the PC1-LD50 comparison)
#   PART 2 (local):  NJ tree + Mantel tests -- produces the genome-wide NJ
#                     tree (Figure S2B) and the CDS-vs-genome-wide
#                     concordance statistic (Mantel r = 0.57, p < 0.001)
#
# Run Part 1 on Savio first (submit via sbatch, ~9 hours), download its
# output, then run Part 2 locally in RStudio.
# =============================================================================


# =============================================================================
# PART 1 (SAVIO) -- genome-wide PCA
#
# Submit as an sbatch job (not interactively -- this needs 8-24 hours):
#   sbatch 05a_genomewide_pca.sh
# Requires: merged_all_samples.vcf.gz, ld50_calculated.csv
# Produces: genomewide_pca_scores.csv, genomewide_pca_result.RData,
#           genomewide_genlight.RData (checkpoint, in case of rerun)
#
# The block below is the R portion of that sbatch script (originally
# genomewide_pca_v2.sh). SBATCH directives for reference:
#   #SBATCH --partition=savio2_bigmem
#   #SBATCH --cpus-per-task=8
#   #SBATCH --time=24:00:00
# =============================================================================

# --- Begin Savio-only R code ---
#
# library(vcfR)
# library(adegenet)
#
# genlight_cache <- "/global/scratch/users/pswhite/analysis/genomewide_genlight.RData"
#
# if (file.exists(genlight_cache)) {
#   cat(format(Sys.time()), "- Loading cached genlight object (skipping VCF re-read)...\n")
#   load(genlight_cache)
# } else {
#   cat(format(Sys.time()), "- Reading genome-wide VCF...\n")
#   vcf <- read.vcfR("/global/scratch/users/pswhite/analysis/merged_all_samples.vcf.gz", verbose = TRUE)
#
#   cat(format(Sys.time()), "- Converting to genlight...\n")
#   genlight_obj <- vcfR2genlight(vcf)
#
#   il_map <- c(I1="IL-1", I2="IL-2", I3="IL-4", I4="IL-7", I5="IL-9", I6="IL-10",
#               I7="IL-14", I8="IL-17", I9="IL-19", I10="IL-24", I11="IL-29", I12="IL-36")
#   bam_id <- sub(".*/(I[0-9]+)\\.bam$", "\\1", indNames(genlight_obj))
#   indNames(genlight_obj) <- unname(il_map[bam_id])
#   stopifnot(!anyNA(indNames(genlight_obj)))
#
#   cat(format(Sys.time()), "- Saving genlight object as checkpoint...\n")
#   save(genlight_obj, file = genlight_cache)
# }
#
# cat(format(Sys.time()), "- Running PCA (glPca, nf=3, parallel, 8 cores)...\n")
# pca_result <- glPca(genlight_obj, nf = 3, parallel = TRUE, n.cores = 8)
# cat(format(Sys.time()), "- PCA complete!\n")
#
# save(genlight_obj, pca_result, file = "/global/scratch/users/pswhite/analysis/genomewide_pca_result.RData")
#
# pca_df <- data.frame(
#   sample = indNames(genlight_obj),
#   PC1 = pca_result$scores[,1],
#   PC2 = pca_result$scores[,2],
#   PC3 = pca_result$scores[,3]
# )
# write.csv(pca_df, "/global/scratch/users/pswhite/analysis/genomewide_pca_scores.csv", row.names=FALSE)
#
# ld50 <- read.csv("/global/scratch/users/pswhite/analysis/ld50_calculated.csv")
# ld50 <- setNames(ld50$Avg_LD50[ld50$Population != "Stock"], ld50$Population[ld50$Population != "Stock"])
# pca_df$ld50 <- ld50[pca_df$sample]
# print(cor.test(pca_df$PC1, pca_df$ld50))
#
# --- End Savio-only R code ---


# =============================================================================
# PART 2 (LOCAL) -- NJ tree + Mantel tests
#
# Run in RStudio after downloading genetic_distances_genomewide_FILTERED.csv
# from Savio (output of the genome-wide distance-matrix pipeline).
# =============================================================================

library(ape)
library(vegan)
library(dplyr)

# ── Load data ────────────────────────────────────────────────────────────────
dist_genomewide <- read.csv("~/Desktop/genetic_distances_genomewide_FILTERED.csv",
                             row.names = 1, check.names = FALSE)
dist_genomewide <- as.matrix(dist_genomewide)

dist_cds <- read.csv("~/Desktop/IL Genomics/genetic_distances.csv",
                      row.names = 1, check.names = FALSE)
dist_cds <- as.matrix(dist_cds)

# Ensure same line order in both matrices
line_order <- rownames(dist_cds)
dist_genomewide <- dist_genomewide[line_order, line_order]

ld50_table <- read.csv("~/Desktop/IL Genomics/ld50_calculated.csv")
ld50_values <- setNames(ld50_table$Avg_LD50, ld50_table$Population)
ld50_values <- ld50_values[names(ld50_values) != "Stock"]
ld50_values <- ld50_values[line_order]

# ── 1. NJ tree from genome-wide, filtered, corrected distances ────────────────
tree_gw <- nj(as.dist(dist_genomewide))

color_ramp <- colorRampPalette(c("steelblue", "gold", "firebrick"))
n_colors   <- 100
colors     <- color_ramp(n_colors)
tip_ld50    <- ld50_values[tree_gw$tip.label]
ld50_scaled <- as.integer((tip_ld50 - min(tip_ld50)) /
                            (max(tip_ld50) - min(tip_ld50)) * (n_colors - 1)) + 1
tip_colors  <- colors[ld50_scaled]

pdf("~/Desktop/Figure_3_NJ_tree_genomewide.pdf", width = 7, height = 8)
par(mar = c(5, 1, 2, 4))
plot(tree_gw, type = "phylogram", tip.color = tip_colors,
     edge.width = 2.5, cex = 1.2, font = 3, direction = "rightwards", no.margin = FALSE)
add.scale.bar(cex = 0.85, col = "gray40")
dev.off()

# ── 2. Mantel test: genome-wide genomic distance vs LD50 distance ─────────────
ld50_dist <- dist(ld50_values, method = "euclidean")
mantel_gw <- mantel(as.dist(dist_genomewide), ld50_dist,
                     method = "pearson", permutations = 9999)
print(mantel_gw)

# ── 3. Mantel test: CDS-only genomic distance vs LD50 distance (original) ─────
mantel_cds <- mantel(as.dist(dist_cds), ld50_dist,
                      method = "pearson", permutations = 9999)
print(mantel_cds)

# ── 4. Does CDS-restriction itself change the structure? ──────────────────────
mantel_cds_vs_gw <- mantel(as.dist(dist_cds), as.dist(dist_genomewide),
                            method = "pearson", permutations = 9999)
print(mantel_cds_vs_gw)

upper_cds <- dist_cds[upper.tri(dist_cds)]
upper_gw  <- dist_genomewide[upper.tri(dist_genomewide)]
cat("Pearson correlation between CDS-only and genome-wide pairwise distances:",
    round(cor(upper_cds, upper_gw), 3), "\n")
