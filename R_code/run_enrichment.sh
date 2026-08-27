#!/bin/bash
# Runs the immune-gene enrichment permutation test end to end.
# Usage: bash run_enrichment.sh
set -e

ANALYSIS=/global/scratch/users/pswhite/analysis
GFF=/global/scratch/users/pswhite/data/plodia_ref/genomic.gff
cd $ANALYSIS

module load r/4.4.0-gcc-11.4.0
module load bio/bedtools2/2.31.0-gcc-11.4.0

echo "STEP 1: writing BED of all analyzed windows..."
Rscript -e '
.libPaths("~/R/library")
options(scipen = 999)   # FIX: prevent scientific notation (1.1e+07) in BED output
setwd("/global/scratch/users/pswhite/analysis")
d <- read.csv("joint_results_FINAL.csv", stringsAsFactors=FALSE)
d$chrom <- sub("_[0-9]+_[0-9]+\\.vcf$", "", d$clean_window)
d$start <- as.integer(sub(".*_([0-9]+)_[0-9]+\\.vcf$", "\\1", d$clean_window))
d$end   <- as.integer(sub(".*_[0-9]+_([0-9]+)\\.vcf$", "\\1", d$clean_window))
write.table(d[,c("chrom","start","end","clean_window")], "all_analyzed_windows.bed",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
cat("wrote", nrow(d), "windows\n")
'

echo "  sanity check -- first 3 BED lines (coordinates must be plain integers):"
head -3 all_analyzed_windows.bed
if grep -q "e+" all_analyzed_windows.bed; then
  echo "ERROR: scientific notation still present in BED file. Stopping."
  exit 1
fi

echo "STEP 2: intersecting all windows with annotation (this takes a minute)..."
bedtools intersect -a all_analyzed_windows.bed -b $GFF -wa -wb \
  | awk -F'\t' '$7=="gene"' \
  | awk -F'\t' '{split($13,a,"description="); split(a[2],b,";"); print $4"\t"b[1]}' \
  > window_gene_table.txt

N_GENES=$(wc -l < window_gene_table.txt)
N_WIN=$(cut -f1 window_gene_table.txt | sort -u | wc -l)
echo "  gene records mapped: $N_GENES  across  $N_WIN  windows"

if [ "$N_WIN" -lt 2000 ]; then
  echo "ERROR: expected ~2164 windows with genes, got $N_WIN. Stopping."
  exit 1
fi

echo "STEP 3: permutation test..."
Rscript $ANALYSIS/enrichment_test.R
