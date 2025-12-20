# 08 - FA VS URBANITY MODEL -----------------------------------------------
#
# Purpose:
#   - Test whether fluctuating asymmetry (FA) differs between
#     urban and rural environments
#
# Model:
#   FA ~ urbanity + (1 | site)
#
# Inputs:
#   results/05_fa_analysis/fa_table_<species>.csv
#
# Outputs:
#   results/08_fa_models/
#     - model_summary_<species>.txt
#     - model_coefficients_<species>.csv
#     - fa_by_urbanity_<species>.png
#     - diagnostics_<species>.png
# ------------------------------------------------------------------------

# 8.1 - Libraries ---------------------------------------------------------

if (!requireNamespace("lme4", quietly = TRUE)) install.packages("lme4")
if (!requireNamespace("lmerTest", quietly = TRUE)) install.packages("lmerTest")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(lme4)
library(lmerTest)
library(ggplot2)

# 8.2 - Settings ----------------------------------------------------------

fa_dir  <- file.path("results", "05_fa_analysis")
out_dir <- file.path("results", "08_fa_models")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Optional: log-transform FA (recommended = TRUE)
log_FA <- TRUE

# 8.3 - Detect species automatically -------------------------------------

fa_files <- list.files(
  fa_dir,
  pattern = "^fa_table_.*\\.csv$",
  full.names = TRUE
)

if (length(fa_files) == 0) {
  stop("No FA tables found in ", fa_dir)
}

species <- gsub("^fa_table_|\\.csv$", "", basename(fa_files))

cat("Species detected:\n")
for (sp in species) cat("  - ", sp, "\n", sep = "")

# 8.4 - Loop over species -------------------------------------------------

for (sp in species) {
  
  cat("\nAnalyzing Bombus", sp, "\n")
  
  fa_file <- file.path(fa_dir, paste0("fa_table_", sp, ".csv"))
  fa <- read.csv(fa_file, stringsAsFactors = FALSE)
  
  # Factors
  fa$urbanity <- factor(fa$urbanity, levels = c("Rural", "Urban"))
  fa$site     <- factor(fa$site)
  
  cat("  Sample size by environment:\n")
  print(table(fa$urbanity))
  
  # Response variable
  fa$response <- if (log_FA) log(fa$FA_unsigned) else fa$FA_unsigned
  response_label <- if (log_FA) "log(FA)" else "FA"
  
  # Model
  m <- lmer(
    response ~ urbanity + (1 | site),
    data = fa,
    REML = TRUE
  )
  
  sm <- summary(m)
  
  # Save model summary
  sink(file.path(out_dir, paste0("model_summary_", sp, ".txt")))
  cat("Bombus", sp, "\n\n")
  cat("Response:", response_label, "\n\n")
  print(sm)
  sink()
  
  # Coefficients
  coef_df <- as.data.frame(sm$coefficients)
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL
  
  write.csv(
    coef_df,
    file.path(out_dir, paste0("model_coefficients_", sp, ".csv")),
    row.names = FALSE
  )
  
  # Plot FA by urbanity
  p <- ggplot(fa, aes(x = urbanity, y = response)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 1.8) +
    labs(
      title = paste("Fluctuating asymmetry – Bombus", sp),
      x = "Environment",
      y = response_label
    ) +
    theme_classic(base_size = 14)
  
  ggsave(
    filename = file.path(out_dir, paste0("fa_by_urbanity_", sp, ".png")),
    plot = p,
    width = 6,
    height = 4,
    dpi = 300
  )
  
  # Diagnostics
  png(
    file.path(out_dir, paste0("diagnostics_", sp, ".png")),
    width = 1600, height = 1600, res = 300
  )
  par(mfrow = c(1, 2))
  plot(fitted(m), resid(m),
       xlab = "Fitted values", ylab = "Residuals",
       main = "Residuals vs fitted")
  abline(h = 0, col = "red")
  qqnorm(resid(m))
  qqline(resid(m), col = "red")
  dev.off()
  
  cat("  Model, plot, and diagnostics saved for Bombus", sp, "\n")
}

cat("\nOutputs written to:\n  ", out_dir, "\n", sep = "")
