#!/bin/bash
#SBATCH --job-name=bin_ld_decay
#SBATCH --account=fc_sirmodel
#SBATCH --partition=savio2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=/global/scratch/users/pswhite/analysis/bin_ld_decay_%j.log

ANALYSIS=/global/scratch/users/pswhite/analysis
INPUT=$ANALYSIS/ld_decay.geno.ld
OUTPUT=$ANALYSIS/ld_decay_binned.csv

# Single streaming pass through the 3.5GB pairwise LD file using awk --
# never loads the full file into memory, just accumulates running sums
# per distance bin. Bins are 50kb wide, out to 2Mb (matching the
# --ld-window-bp cap used when the file was generated).
#
# Input columns (tab-separated, from vcftools --geno-r2):
#   CHR  POS1  POS2  N_INDV  R^2
#
# Skips pairs where R^2 is nan/-nan (vcftools outputs this when a site
# pair can't be evaluated, e.g. monomorphic after filtering).

echo "Binning pairwise LD by distance..."

awk -F'\t' '
    NR == 1 { next }  # skip header row
    {
        r2 = $5
        if (r2 == "nan" || r2 == "-nan") next
        dist = $3 - $2
        if (dist < 0) dist = -dist
        bin = int(dist / 50000)
        sum[bin] += r2
        n[bin] += 1
    }
    END {
        print "bin_start_bp,bin_end_bp,mean_r2,n_pairs"
        for (b = 0; b <= 40; b++) {
            if (n[b] > 0) {
                printf "%d,%d,%.6f,%d\n", b*50000, (b+1)*50000, sum[b]/n[b], n[b]
            }
        }
    }
' "$INPUT" > "$OUTPUT"

echo "Done. Output written to $OUTPUT"
echo "Preview:"
cat "$OUTPUT"
