library(dplyr)
library(mclust)   # for proper BIC-based model selection with small n

pca <- read.csv("~/Desktop/IL Genomics/cds_pca_scores_FRESH.csv")
X <- as.matrix(pca[, c("PC1","PC2","PC3")])
rownames(X) <- pca$sample

# ── WSS across k = 1 to 5 (same diagnostic style as Table S7/S8's LD50 clustering) ──
wss <- sapply(1:5, function(k) {
  if (k == 1) return(sum(scale(X, scale=FALSE)^2))
  set.seed(42)
  kmeans(X, centers=k, nstart=25)$tot.withinss
})
cat("WSS by k:\n")
print(round(wss, 2))
cat("\n% improvement at each step:\n")
print(round((wss[-5] - wss[-1]) / wss[-5] * 100, 1))

# ── Gaussian mixture model selection via BIC (mclust) -- more rigorous for n=12 ──
mc <- Mclust(X, G=1:5)
cat("\nBest model by BIC:\n")
print(summary(mc))
cat("\nOptimal number of clusters (BIC-selected):", mc$G, "\n")

# ── Silhouette check as a third, independent diagnostic ──
library(cluster)
sil_widths <- sapply(2:5, function(k) {
  set.seed(42)
  km <- kmeans(X, centers=k, nstart=25)
  mean(silhouette(km$cluster, dist(X))[,3])
})
cat("\nAverage silhouette width, k=2 to 5:\n")
print(round(sil_widths, 3))

cat("Per-sample cluster assignment (Mclust, G=5):\n")
print(data.frame(sample = rownames(X), cluster = mc$classification))

cat("PC1 alone, sorted:\n")
print(sort(pca$PC1))
cat("\nHartigan's dip test for unimodality (PC1 only):\n")
library(diptest)
print(dip.test(pca$PC1))
