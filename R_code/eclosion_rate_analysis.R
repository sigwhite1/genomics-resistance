# =============================================================================
# eclosion_rate_analysis.R
#
# Calculates per-population eclosion rate (proportion of pupated individuals
# that successfully eclosed) from raw pupal-weight tracking data, and tests
# whether it differs among inbred lines and correlates with resistance.
# Addresses Reviewer 2's comment on eclosion rate as a life-history trait.
#
# Input:  Pupal_Weights_Data_II.xlsx (one sheet per population; per-individual
#         pupation date and eclosion outcome, recorded either as a "Date of
#         Eclosion" column or a free-text "eclosed" annotation depending on
#         population -- see notes below)
#         ld50_calculated.csv
#
# Output: Table S4 eclosion columns (Pupated, Eclosed, Eclosion rate)
#         Results-section statistics: chi-square test across lines,
#         Pearson correlation with LD50
#
# NOTE: IL-1's second assay round (46 additional pupae, tracked separately
# from the first round in the same sheet) is excluded, because eclosion
# outcome was not recorded for that round at all (no "eclosed" annotation,
# no Date of Eclosion column, and no columns beyond weight measurements).
# Including it would incorrectly count those individuals as non-eclosed.
# =============================================================================

library(readxl)
library(dplyr)

lines <- c("IL-1","IL-2","IL-4","IL-7","IL-9","IL-10",
           "IL-14","IL-17","IL-19","IL-24","IL-29","IL-36","Stock")

xl <- "~/Desktop/IL Genomics/Pupal_Weights_Data_II.xlsx"

results <- list()

for (s in lines) {
  df <- read_excel(xl, sheet = s)
  df <- df[!is.na(df$Sex), ]

  if (s == "IL-1") {
    # Exclude untracked second round (Sex recorded as "Females"/"Males",
    # plural, vs. "Female"/"Male" for the first, tracked round)
    df <- df[!grepl("s$", df$Sex), ]
  }

  pup_col <- grep("Date of Pupation", names(df), value = TRUE)[1]
  pupated <- df[!is.na(df[[pup_col]]), ]

  ecl_col <- grep("Eclosion", names(df), value = TRUE)
  if (length(ecl_col) > 0) {
    eclosed <- sum(!is.na(pupated[[ecl_col[1]]]))
  } else {
    eclosed <- sum(apply(pupated, 1, function(row)
      any(grepl("eclos", row, ignore.case = TRUE))))
  }

  results[[s]] <- data.frame(Population = s, Pupated = nrow(pupated), Eclosed = eclosed)
}

eclosion_table <- bind_rows(results) %>%
  mutate(Eclosion_Rate = Eclosed / Pupated)

print(eclosion_table)

# ── Chi-square test of independence across the 12 inbred lines ───────────────
inbred <- eclosion_table %>% filter(Population != "Stock")
chisq_tab <- as.matrix(inbred[, c("Eclosed", "Pupated")] %>%
                          mutate(Not_Eclosed = Pupated - Eclosed) %>%
                          select(Eclosed, Not_Eclosed))
rownames(chisq_tab) <- inbred$Population
print(chisq.test(chisq_tab))

# ── Correlation with LD50 ─────────────────────────────────────────────────────
ld50_table <- read.csv("~/Desktop/IL Genomics/ld50_calculated.csv")
merged <- merge(inbred, ld50_table, by.x = "Population", by.y = "Population")
print(cor.test(merged$Eclosion_Rate, merged$Avg_LD50))
