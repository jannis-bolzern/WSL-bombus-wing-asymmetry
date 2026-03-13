# 08 - FA VS URBANITY MODEL -----------------------------------------------
#
# Purpose:
#   - Test whether fluctuating asymmetry (FA) differs between
#     urban and rural environments
#
# Model:
#   response ~ urbanity + (1 | site)
#
# Inputs:
#   results/05_extract_fa/fa_table_<species>.csv
#
# Outputs:
#   results/08_fa_models/
#     - model_summary_<species>.txt
#     - model_coefficients_<species>.csv
#     - model_effect_<species>.csv              (NEW: effect + CI)
#     - fa_by_urbanity_<species>.png            (improved plot)
#     - diagnostics_<species>.png

# 8.1 - Libraries ---------------------------------------------------------

if (!requireNamespace("lme4", quietly = TRUE)) install.packages("lme4")
if (!requireNamespace("lmerTest", quietly = TRUE)) install.packages("lmerTest")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(lme4)
library(lmerTest)
library(ggplot2)

# 8.2 - Settings ----------------------------------------------------------

fa_dir  <- file.path("results", "05_extract_fa")
out_dir <- file.path("results", "08_fa_models")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 8.3 - Detect species automatically --------------------------------------

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
  
  cat("  Sample size by environment (individuals):\n")
  print(table(fa$urbanity))
  
  # ---- Response variable ----
  fa$response <- fa$FA_unsigned
  response_label <- "Unsigned FA (Procrustes units)"
  
  # ---- Model ----
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
  
  # Coefficients table
  coef_df <- as.data.frame(sm$coefficients)
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL
  
  write.csv(
    coef_df,
    file.path(out_dir, paste0("model_coefficients_", sp, ".csv")),
    row.names = FALSE
  )
  
  # Effect size + CI (Urban - Rural)
  beta <- fixef(m)["urbanityUrban"]
  ci <- suppressMessages(confint(m, parm = "urbanityUrban", method = "Wald"))
  pval <- coef(summary(m))["urbanityUrban", "Pr(>|t|)"]
  
  eff <- data.frame(
    species = sp,
    response = response_label,
    beta_urban_minus_rural = as.numeric(beta),
    ci_low = as.numeric(ci[1]),
    ci_high = as.numeric(ci[2]),
    p_value = as.numeric(pval),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    eff,
    file.path(out_dir, paste0("model_effect_", sp, ".csv")),
    row.names = FALSE
  )
  
  # --- Boxplot ---
  
  n_tab <- table(fa$urbanity)
  x_labs <- c(
    Rural = paste0("Rural\n(n=", n_tab["Rural"], ")"),
    Urban = paste0("Urban\n(n=", n_tab["Urban"], ")")
  )
  
  ann <- sprintf("p = %.3g", as.numeric(pval))
  
  y_max <- max(fa$response, na.rm = TRUE)
  y_rng <- diff(range(fa$response, na.rm = TRUE))
  y_br  <- y_max + 0.06 * y_rng
  y_txt <- y_max + 0.10 * y_rng
  
  p <- ggplot(fa, aes(x = urbanity, y = response)) +
    geom_boxplot(
      width = 0.55,
      linewidth = 0.7,
      outlier.size = 1.3,
      outlier.alpha = 0.45
    ) +
    # bracket (no warnings)
    annotate("segment", x = 1, xend = 2, y = y_br, yend = y_br, linewidth = 0.6) +
    annotate("segment", x = 1, xend = 1, y = y_br, yend = y_br - 0.02*y_rng, linewidth = 0.6) +
    annotate("segment", x = 2, xend = 2, y = y_br, yend = y_br - 0.02*y_rng, linewidth = 0.6) +
    # effect text
    annotate("text", x = 1.5, y = y_txt, label = ann, size = 3.6, vjust = 0) +
    scale_x_discrete(labels = x_labs) +
    labs(
      title = paste("FA vs urbanity – Bombus", sp),
      x = NULL,
      y = response_label
    ) +
    coord_cartesian(ylim = c(min(fa$response, na.rm = TRUE), y_max + 0.18*y_rng)) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(lineheight = 0.95)
    )
  
  ggsave(
    filename = file.path(out_dir, paste0("fa_by_urbanity_", sp, ".png")),
    plot = p,
    width = 5.8,
    height = 4.3,
    dpi = 300
  )
  
  # ---- Diagnostics ----
  png(
    file.path(out_dir, paste0("diagnostics_", sp, ".png")),
    width = 1600, height = 1600, res = 300
  )
  par(mfrow = c(2, 2))
  plot(fitted(m), resid(m),
       xlab = "Fitted values", ylab = "Residuals",
       main = "Residuals vs fitted")
  abline(h = 0, col = "red")
  qqnorm(resid(m)); qqline(resid(m), col = "red")
  hist(resid(m), main = "Residual histogram", xlab = "Residuals")
  dev.off()
  
  cat("  Model, effect table, plot, and diagnostics saved for Bombus", sp, "\n")
}

cat("\nOutputs written to:\n  ", out_dir, "\n", sep = "")
