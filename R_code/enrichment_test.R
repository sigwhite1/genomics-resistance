# =============================================================================
# enrichment_test.R
#
# Permutation test: are the peak convergent windows enriched for immune-
# related genes relative to random windows from the same analyzed set?
#
# NULL: draw N random windows (N = number of peak windows) from the 2164
# analyzed windows, 10,000 times, and count immune-keyword matches. Compare
# the observed count against that distribution.
#
# CAVEAT ON KEYWORDS: this list was written after the peak-window gene list
# had already been read, which is a contamination risk. It is deliberately
# comprehensive -- drawn from canonical insect immunity pathways and
# including many terms absent from the observed list -- so that any
# unconscious tuning toward observed genes works against significance
# rather than for it. Strict and broad variants are both reported.
# =============================================================================

.libPaths("~/R/library")
setwd("/global/scratch/users/pswhite/analysis")
set.seed(42)

tab <- read.delim("window_gene_table.txt", header=FALSE,
                  col.names=c("window","description"), stringsAsFactors=FALSE)
tab$description[is.na(tab$description)] <- ""

peaks <- read.csv("convergent_regions_summary.csv", stringsAsFactors=FALSE)$best_window
cat("Peak windows:", length(peaks), "\n")
cat("Total windows with genes:", length(unique(tab$window)), "\n")
cat("Total gene records:", nrow(tab), "\n\n")

# --- Pre-specified keyword sets -------------------------------------------
strict <- paste(c(
  "toll","spatzle","spaetzle","cactus","pelle","myd88",
  "\\bimd\\b","relish","caspar","dredd",
  "\\bjak\\b","\\bstat\\b","hopscotch","domeless","unpaired",
  "argonaute","dicer","drosha","loquacious","r2d2",
  "phenoloxidase","melaniz","serpin","clip-domain",
  "lysozyme","defensin","cecropin","attacin","gloverin","moricin",
  "lebocin","diptericin","drosomycin","metchnikowin","gambicin",
  "peptidoglycan","\\bpgrp\\b","gram-negative bacteria-binding","beta-1,3-glucan",
  "thioester-containing","complement",
  "hemocyte","phagocyt","encapsulat","nodulat",
  "nimrod","eater","draper","croquemort","scavenger receptor",
  "dscam","down syndrome cell adhesion",
  "cyclic gmp-amp","sting","cgas",
  "immune","immunity","antimicrobial","antiviral","antibacterial","defense response",
  "transferrin","ferritin","nitric oxide synthase","hemolin",
  "gamma-interferon-inducible","thiol reductase"
), collapse="|")

broad_extra <- paste(c(
  "proteasome","ubiquitin","deubiquitinase",
  "autophag","\\batg[0-9]","light chain 3","unc-51","\\bulk\\b","beclin",
  "dual oxidase","\\bduox\\b","nadph oxidase","superoxide dismutase",
  "catalase","peroxidase","peroxiredoxin","glutathione s-transferase","thioredoxin",
  "heat shock","chaperone",
  "caspase","apoptos",
  "cytochrome p450","carboxylesterase",
  "serine protease","trypsin","chymotrypsin","serine-type endopeptidase"
), collapse="|")
broad <- paste(strict, broad_extra, sep="|")

run_test <- function(pattern, label) {
  tab$hit <- grepl(pattern, tab$description, ignore.case=TRUE)
  per_window <- tapply(tab$hit, tab$window, sum)
  all_windows <- names(per_window)

  obs <- sum(per_window[names(per_window) %in% peaks], na.rm=TRUE)
  n_draw <- sum(peaks %in% all_windows)

  null_counts <- replicate(10000, sum(per_window[sample(all_windows, n_draw)]))
  p <- (sum(null_counts >= obs) + 1) / (10001)

  cat("=== ", label, " ===\n", sep="")
  cat("Observed immune-keyword genes in peak windows:", obs, "\n")
  cat("Null mean:", round(mean(null_counts),2),
      " sd:", round(sd(null_counts),2),
      " 95th pctile:", quantile(null_counts, 0.95), "\n")
  cat("Permutation p-value (one-sided):", signif(p,4), "\n")
  cat("Genome-wide keyword rate:", round(100*mean(tab$hit),2), "% of genes\n\n")

  cat("Matching genes in peak windows:\n")
  m <- tab[tab$window %in% peaks & tab$hit, c("window","description")]
  if (nrow(m)) print(m, row.names=FALSE) else cat("  (none)\n")
  cat("\n")
}

run_test(strict, "STRICT: canonical immune pathway genes")
run_test(broad,  "BROAD: immune + proteostasis/autophagy/oxidative stress/proteases")
