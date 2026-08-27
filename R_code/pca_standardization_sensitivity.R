# =============================================================================
# pca_standardization_sensitivity.R
#
# Tests whether the genomic PC1-LD50 association (the manuscript's headline
# result) is sensitive to the choice of PCA standardization. glPca's default
# does not standardize SNPs by allele frequency before PCA; this script
# repeats the analysis with standardization enabled (scale=TRUE), the
# convention used by EIGENSTRAT/smartpca, restricted to the same marker set
# used by an independent check via SNPRelate (which standardizes by default).
#
# Addresses Reviewer 1's comment on sample-size/methodological sensitivity of
# the headline result. Reported in Results/Discussion as r = -0.58, p = 0.049
# (unstandardized: r = -0.65, p = 0.022).
#
# Run on Savio HPC (requires the CDS-restricted VCF and a completed SNPRelate
# run to identify the marker subset).
# =============================================================================

library(SNPRelate)
library(vcfR)
library(adegenet)

setwd("/global/scratch/users/pswhite/data/inbred_only/bcf/")

# ── Step 1: get the exact marker set SNPRelate retained ────────────────────
# (SNPRelate excludes monomorphic/zero-MAF/all-missing SNPs by default;
# restricting glPca to the same set isolates standardization as the only
# remaining difference between the two methods.)
genofile <- snpgdsOpen("cds_verify.gds")
pca_check <- snpgdsPCA(genofile, autosome.only = FALSE)
used_snp_id <- pca_check$snp.id
all_snp_id  <- read.gdsn(index.gdsn(genofile, "snp.id"))
kept_chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
kept_pos <- read.gdsn(index.gdsn(genofile, "snp.position"))
keep_idx <- match(used_snp_id, all_snp_id)
used_chr <- kept_chr[keep_idx]
used_pos <- kept_pos[keep_idx]
snpgdsClose(genofile)

# ── Step 2: subset the genlight object to the same markers ─────────────────
vcf <- read.vcfR("merged_inbred_intersect_cds.vcf.gz", verbose = FALSE)
genlight_obj <- vcfR2genlight(vcf)
gl_chr <- vcf@fix[,"CHROM"]
gl_pos <- as.integer(vcf@fix[,"POS"])
match_key <- paste(used_chr, used_pos, sep = "_")
vcf_key   <- paste(gl_chr, gl_pos, sep = "_")
keep_loci <- vcf_key %in% match_key
genlight_subset <- genlight_obj[, keep_loci]

il_map <- c(I1="IL-1", I2="IL-2", I3="IL-4", I4="IL-7", I5="IL-9", I6="IL-10",
            I7="IL-14", I8="IL-17", I9="IL-19", I10="IL-24", I11="IL-29", I12="IL-36")
bam_id <- sub(".*/(I[0-9]+)\\.bam$", "\\1", indNames(genlight_subset))
indNames(genlight_subset) <- unname(il_map[bam_id])

ld50 <- read.csv("/global/scratch/users/pswhite/analysis/ld50_calculated.csv")
ld50 <- setNames(ld50$Avg_LD50[ld50$Population != "Stock"],
                  ld50$Population[ld50$Population != "Stock"])

# ── Step 3: sanity check -- unscaled glPca on the matched subset ───────────
# (should reproduce the same magnitude/significance as the full CDS analysis,
# confirming the marker-subsetting itself isn't driving any difference)
cat("=== glPca(scale=FALSE, default) on SNPRelate-matched marker subset ===\n")
pca_unscaled <- glPca(genlight_subset, nf = 3, scale = FALSE)
pc1_u <- setNames(pca_unscaled$scores[,1], indNames(genlight_subset))
print(cor.test(pc1_u, ld50[names(pc1_u)]))

# ── Step 4: the actual sensitivity check -- standardized glPca ─────────────
cat("\n=== glPca(scale=TRUE) on SNPRelate-matched marker subset ===\n")
pca_scaled <- glPca(genlight_subset, nf = 3, scale = TRUE)
pc1_s <- setNames(pca_scaled$scores[,1], indNames(genlight_subset))
print(cor.test(pc1_s, ld50[names(pc1_s)]))

# ── Step 5 (context only, not part of this sensitivity check itself) ──────
# An independent third method (SNPRelate's own PCA, which standardizes by
# allele frequency by default) on the same marker subset:
cat("\n=== snpgdsPCA (SNPRelate, same marker subset) -- for reference ===\n")
genofile <- snpgdsOpen("cds_verify.gds")
pca_snprelate <- snpgdsPCA(genofile, autosome.only = FALSE)
bam_id2 <- sub(".*/(I[0-9]+)\\.bam$", "\\1", pca_snprelate$sample.id)
pc1_snprelate <- setNames(pca_snprelate$eigenvect[,1], unname(il_map[bam_id2]))
print(cor.test(pc1_snprelate, ld50[names(pc1_snprelate)]))
snpgdsClose(genofile)
