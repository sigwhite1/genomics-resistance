#!/bin/bash
#SBATCH --job-name=genomewide_pca_v2
#SBATCH --account=fc_sirmodel
#SBATCH --partition=savio2_bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=/global/scratch/users/pswhite/analysis/genomewide_pca_v2_%j.log

module load r/4.4.0-gcc-11.4.0

Rscript - <<'EOF'
library(vcfR)
library(adegenet)

genlight_cache <- "/global/scratch/users/pswhite/analysis/genomewide_genlight.RData"

if (file.exists(genlight_cache)) {
  cat(format(Sys.time()), "- Loading cached genlight object (skipping VCF re-read)...\n")
  load(genlight_cache)
} else {
  cat(format(Sys.time()), "- Reading genome-wide VCF...\n")
  vcf <- read.vcfR("/global/scratch/users/pswhite/analysis/merged_all_samples.vcf.gz", verbose = TRUE)

  cat(format(Sys.time()), "- Converting to genlight...\n")
  genlight_obj <- vcfR2genlight(vcf)

  il_map <- c(I1="IL-1", I2="IL-2", I3="IL-4", I4="IL-7", I5="IL-9", I6="IL-10",
              I7="IL-14", I8="IL-17", I9="IL-19", I10="IL-24", I11="IL-29", I12="IL-36")
  bam_id <- sub(".*/(I[0-9]+)\\.bam$", "\\1", indNames(genlight_obj))
  indNames(genlight_obj) <- unname(il_map[bam_id])
  stopifnot(!anyNA(indNames(genlight_obj)))

  cat(format(Sys.time()), "- Saving genlight object as checkpoint...\n")
  save(genlight_obj, file = genlight_cache)
}

cat(format(Sys.time()), "- Running PCA (glPca, nf=3, parallel, 8 cores)...\n")
pca_result <- glPca(genlight_obj, nf = 3, parallel = TRUE, n.cores = 8)
cat(format(Sys.time()), "- PCA complete!\n")

save(genlight_obj, pca_result, file = "/global/scratch/users/pswhite/analysis/genomewide_pca_result.RData")

cat("\nVariance explained (%):\n")
print(round(pca_result$eig / sum(pca_result$eig) * 100, 1)[1:5])

pca_df <- data.frame(
  sample = indNames(genlight_obj),
  PC1 = pca_result$scores[,1],
  PC2 = pca_result$scores[,2],
  PC3 = pca_result$scores[,3]
)
write.csv(pca_df, "/global/scratch/users/pswhite/analysis/genomewide_pca_scores.csv", row.names=FALSE)

cat("\nGenomic PC1 - LD50 correlation:\n")
ld50 <- read.csv("/global/scratch/users/pswhite/analysis/ld50_calculated.csv")
ld50 <- setNames(ld50$Avg_LD50[ld50$Population != "Stock"], ld50$Population[ld50$Population != "Stock"])
pca_df$ld50 <- ld50[pca_df$sample]
ct <- cor.test(pca_df$PC1, pca_df$ld50)
print(ct)

cat(format(Sys.time()), "- Done!\n")
EOF
