# =============================================================================
# 05_ROH_figure.R
#
# Purpose: Supplementary figure showing runs of homozygosity (ROH) per
#          inbred line vs stock population, coloured by LD50 for inbred lines.
#
# Inputs:  none (values hardcoded from bcftools roh output on Savio)
#
# Outputs:
#   - FigureS3_ROH.pdf
#   - Console: ROH vs LD50 correlation (Pearson's r and p-value)
#
# NOTE: Both inbred lines and stock "samples" are pooled (2-3 individuals
#       per DNA extraction). n values in figure labels reflect pooled samples,
#       not individual insects.
# =============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)

out_dir <- "~/Desktop/IL Genomics/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

POP_ORDER <- c("IL-1","IL-2","IL-4","IL-7","IL-9","IL-10",
               "IL-14","IL-17","IL-19","IL-24","IL-29","IL-36")

genome_mb <- 291.3

# Inbred line ROH data (each line = 1 pooled sample of 2-3 individuals)
roh_inbred <- data.frame(
  Population = POP_ORDER,
  ROH_Mb     = c(71.1, 82.2, 42.0, 89.8, 61.7, 79.3,
                 73.6, 93.5, 82.3, 82.3, 64.7, 90.9),
  type       = "Inbred line"
)

ld50 <- read.csv("~/Desktop/IL Genomics/ld50_calculated.csv")
ld50 <- ld50[ld50$Population != "Stock", ]
roh_inbred <- merge(roh_inbred, ld50[, c("Population","Avg_LD50")], by="Population")
roh_inbred$Population <- factor(roh_inbred$Population, levels=POP_ORDER)

# Stock ROH data (84 pooled samples, each representing 2-3 individuals)
stock_roh_mb <- c(15.1,9.5,9.8,17.2,13.0,12.9,14.7,31.8,3.7,10.6,5.8,6.8,
                  38.0,12.4,8.3,22.0,9.9,12.5,21.1,22.7,8.3,16.5,8.8,15.3,
                  28.8,18.4,28.8,29.7,19.0,16.1,10.2,32.8,13.0,26.9,28.6,
                  15.2,8.8,32.0,21.6,39.5,25.5,24.7,23.1,25.3,27.4,55.6,
                  13.0,28.3,12.7,33.3,25.7,22.6,28.8,29.0,24.2,35.3,49.5,
                  14.0,16.1,7.8,21.6,25.4,12.1,18.8,14.8,20.0,21.8,27.7,
                  17.8,32.8,21.6,17.0,25.6,26.9,22.0,17.7,26.9,8.3,77.4)

stock_pct <- stock_roh_mb / genome_mb * 100

# Compute % genome in ROH for inbred lines
roh_inbred$pct_genome <- roh_inbred$ROH_Mb / genome_mb * 100
inbred_mean_pct <- mean(roh_inbred$pct_genome)

# =============================================================================
# CORRELATION: ROH vs LD50
# =============================================================================

cor_roh <- cor.test(roh_inbred$pct_genome, roh_inbred$Avg_LD50)
r_roh   <- round(cor_roh$estimate, 2)
p_roh   <- round(cor_roh$p.value, 3)

cat("== ROH vs LD50 correlation ==\n")
cat("r =", r_roh, "  p =", p_roh, "\n")
print(cor_roh)

# =============================================================================
# PANEL A: inbred lines bar chart coloured by LD50
# =============================================================================

p_inbred <- ggplot(roh_inbred, aes(x=Population, y=pct_genome, fill=Avg_LD50)) +
  geom_col(alpha=0.9, width=0.7, color="black", linewidth=0.4) +
  geom_hline(yintercept=inbred_mean_pct, linetype="dashed", color="gray40") +
  annotate("text", x=0.6, y=inbred_mean_pct+0.8,
           label=paste0("mean = ", round(inbred_mean_pct,1), "%"),
           hjust=0, size=3.5, color="gray40") +
  scale_fill_gradientn(colors=c("steelblue","white","firebrick"), name="LD50") +
  scale_y_continuous(limits=c(0,40), breaks=seq(0,40,10), expand=c(0,0)) +
  labs(title="A) Inbred lines", x="Population", y="% genome in ROH") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        panel.border=element_rect(color="black", fill=NA, linewidth=0.5))

# =============================================================================
# PANEL B: stock pooled samples vs inbred line pooled samples
# =============================================================================

compare_df <- data.frame(
  Group   = c(rep("Stock\n(n=84)", length(stock_pct)),
              rep("Inbred lines\n(n=12)", nrow(roh_inbred))),
  pct_ROH = c(stock_pct, roh_inbred$pct_genome)
)
compare_df$Group <- factor(compare_df$Group,
                           levels=c("Stock\n(n=84)", "Inbred lines\n(n=12)"))

p_compare <- ggplot(compare_df, aes(x=Group, y=pct_ROH, fill=Group)) +
  geom_boxplot(alpha=0.8, outlier.shape=NA) +
  geom_jitter(width=0.15, size=1.5, alpha=0.5, color="grey30") +
  scale_fill_manual(values=c("Stock\n(n=84)"="steelblue",
                             "Inbred lines\n(n=12)"="firebrick")) +
  scale_y_continuous(limits=c(0,40), breaks=seq(0,40,10), expand=c(0,0)) +
  labs(title="B) Inbred lines vs. stock population",
       x=NULL, y="% genome in ROH") +
  theme_classic() +
  theme(legend.position="none",
        panel.border=element_rect(color="black", fill=NA, linewidth=0.5))

# =============================================================================
# COMBINE AND SAVE
# =============================================================================

fig_roh <- p_inbred + p_compare + plot_layout(widths=c(2,1))

ggsave(file.path(out_dir, "FigureS3_ROH.pdf"), fig_roh, width=10, height=5)
cat("FigureS3_ROH.pdf saved.\n")
cat("Mean % genome in ROH - Stock:", round(mean(stock_pct),1), "%\n")
cat("Mean % genome in ROH - Inbred:", round(inbred_mean_pct,1), "%\n")
