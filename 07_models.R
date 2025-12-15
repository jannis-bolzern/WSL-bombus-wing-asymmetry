# ============================================================
# 06 — STATISTICAL MODELING OF FA
# ============================================================
# Purpose:
#   - Test environmental effects on FA
#   - Urban vs rural
#   - Species effects
#   - Optional interactions
#
# Input:
#   results/05_FA/FA_per_specimen.csv
#
# Output:
#   results/06_models/
#     - model_summaries.txt
#     - FA_by_urbanity.png
#     - FA_distribution.png
# ============================================================

library(lme4)
library(lmerTest)
library(ggplot2)

# ---------------------------
# Paths
# ---------------------------
fa_file <- file.path("results", "05_FA", "FA_per_specimen.csv")
out_dir <- file.path("results", "06_models")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# Load FA data
# ---------------------------
fa <- read.csv(fa_file, stringsAsFactors = FALSE)

fa$urbanity <- factor(fa$urbanity, levels = c("Rural","Urban"))
fa$species  <- factor(fa$species)
fa$site     <- factor(fa$site)

# ---------------------------
# Quick diagnostics
# ---------------------------
png(file.path(out_dir, "FA_distribution.png"), 800, 600)
hist(fa$FA_unsigned, breaks=20,
     main="Distribution of FA",
     xlab="Unsigned FA")
dev.off()

png(file.path(out_dir, "FA_by_urbanity.png"), 800, 600)
boxplot(FA_unsigned ~ urbanity, data=fa,
        ylab="Unsigned FA",
        main="Fluctuating asymmetry by urbanity")
dev.off()

# ---------------------------
# Models
# ---------------------------

# 1. Simple linear model
m1 <- lm(FA_unsigned ~ urbanity, data = fa)

# 2. Add species
m2 <- lm(FA_unsigned ~ urbanity + species, data = fa)

# 3. Interaction model
m3 <- lm(FA_unsigned ~ urbanity * species, data = fa)

# 4. Mixed model (site as random effect)
m4 <- lmer(FA_unsigned ~ urbanity + species + (1|site), data = fa)

# ---------------------------
# Save summaries
# ---------------------------
sink(file.path(out_dir, "model_summaries.txt"))

cat("Model 1: FA ~ urbanity\n")
print(summary(m1))

cat("\nModel 2: FA ~ urbanity + species\n")
print(summary(m2))

cat("\nModel 3: FA ~ urbanity * species\n")
print(summary(m3))

cat("\nModel 4: FA ~ urbanity + species + (1|site)\n")
print(summary(m4))

sink()

