#!/bin/bash
#SBATCH --job-name=ld_decay
#SBATCH --account=fc_sirmodel
#SBATCH --partition=savio2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=/global/scratch/users/pswhite/analysis/ld_decay_%j.log

module load bio/vcftools/0.1.16-gcc-11.4.0

ANALYSIS=/global/scratch/users/pswhite/analysis
MERGED_VCF=$ANALYSIS/merged_all_samples.vcf.gz

# FIX (2026-07): the original version of this script computed pairwise LD
# across EVERY variant pair within 2Mb, genome-wide, with no thinning.
# This produced >90GB of output in 2+ hours and was cancelled before
# completion. Pairwise LD calculations scale roughly with the square of
# SNP density, so an unthinned genome-wide run is computationally
# infeasible for this purpose.
#
# Fix: thin to one SNP per 2kb before computing LD. This still gives a
# well-spaced, representative sample of SNP pairs across the full range
# of physical distances needed to characterize LD decay, without the
# combinatorial blowup of using every site. The input file remains
# merged_all_samples.vcf.gz (genome-wide, unfiltered) rather than the
# CDS-only file, since this is the same dataset the window-based
# WZA/PCA-LD50 analysis itself uses (see make_windows.sh) -- we want the
# LD-decay estimate to reflect the actual data going into that analysis.
echo "Computing pairwise LD (r^2) up to 2 Mb apart, thinned to 1 SNP per 2kb..."
vcftools --gzvcf $MERGED_VCF \
    --thin 2000 \
    --geno-r2 \
    --ld-window-bp 2000000 \
    --out $ANALYSIS/ld_decay

echo "Done. Output: $ANALYSIS/ld_decay.geno.ld"
echo "Row count: $(wc -l < $ANALYSIS/ld_decay.geno.ld)"
echo "File size: $(du -h $ANALYSIS/ld_decay.geno.ld)"