# 07 - DIGITIZATION ERROR -------------------------------------------------
#
# Measures digitization error using replicate digitizations (series 1 vs 2)
# of the SAME wing (specimen_uid + side).
#
# Outputs:
#   results/07_digitization_error/
#     - digitization_error_per_wing_<SPECIES>.csv
#     - digitization_error_summary_<SPECIES>.csv
# ------------------------------------------------------------------------

# 7.1 Libraries -----------------------------------------------------------

if (!requireNamespace("geomorph", quietly = TRUE)) {
  install.packages("geomorph")
}
library(geomorph)

# 7.2 Settings ------------------------------------------------------------

species <- "lapidarius"   # <- change if needed

data_dir <- "data"
coords_file <- file.path(data_dir, paste0("fw_coords_", species, ".RDS"))
meta_file   <- file.path(data_dir, paste0("fw_metadata_", species, ".csv"))

out_dir <- file.path("results", "07_digitization_error")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 7.3 Load data -----------------------------------------------------------

if (!file.exists(coords_file) || !file.exists(meta_file)) {
  stop("Missing input files for species: ", species)
}

coords_all <- readRDS(coords_file)
meta_all   <- read.csv(meta_file, stringsAsFactors = FALSE)

if (nrow(meta_all) != dim(coords_all)[3]) {
  stop("Metadata / coordinate mismatch for species: ", species)
}

meta_all$side   <- ifelse(meta_all$side %in% c("L","Left"), "Left", "Right")
meta_all$series <- as.integer(meta_all$series)

# wing identity = same physical wing
meta_all$wing_id <- paste(meta_all$specimen_uid, meta_all$side, sep = "__")

# 7.4 Keep only replicated wings -----------------------------------------

tab <- with(meta_all, table(wing_id, series))

if (!all(c("1","2") %in% colnames(tab))) {
  stop("No replicate digitizations found (both series 1 and 2 required).")
}

wing_ids_rep <- rownames(tab)[tab[, "1"] == 1 & tab[, "2"] == 1]

sel <- meta_all$wing_id %in% wing_ids_rep & meta_all$series %in% c(1, 2)
meta   <- meta_all[sel, , drop = FALSE]
coords <- coords_all[,,sel, drop = FALSE]

# sanity check
counts <- table(meta$wing_id, meta$series)
bad <- rownames(counts)[!(counts[, "1"] == 1 & counts[, "2"] == 1)]

if (length(bad) > 0) {
  stop(
    "Some wings do not have exactly one series-1 and one series-2 digitization.\n",
    "Examples:\n",
    paste(head(bad, 10), collapse = ", ")
  )
}

cat("Replicated wings used:", length(unique(meta$wing_id)), "\n")
cat("Digitizations used:", nrow(meta), "\n")

# 7.5 Pairwise Procrustes RMS per wing ------------------------------------

pairwise_procrustes <- function(A, B) {
  
  # center
  A0 <- sweep(A, 2, colMeans(A), "-")
  B0 <- sweep(B, 2, colMeans(B), "-")
  
  # scale to centroid size
  csA <- sqrt(sum(A0^2))
  csB <- sqrt(sum(B0^2))
  
  if (csA > 0) A0 <- A0 / csA
  if (csB > 0) B0 <- B0 / csB
  
  # optimal rotation
  M <- t(B0) %*% A0
  sv <- svd(M)
  R <- sv$u %*% t(sv$v)
  
  if (det(R) < 0) {
    sv$u[, ncol(sv$u)] <- -sv$u[, ncol(sv$u)]
    R <- sv$u %*% t(sv$v)
  }
  
  list(A = A0, B = B0 %*% R)
}

wing_ids <- sort(unique(meta$wing_id))
rms <- numeric(length(wing_ids))

for (i in seq_along(wing_ids)) {
  
  rows <- which(meta$wing_id == wing_ids[i])
  
  i1 <- rows[meta$series[rows] == 1]
  i2 <- rows[meta$series[rows] == 2]
  
  A <- coords[,,i1]
  B <- coords[,,i2]
  
  pr <- pairwise_procrustes(A, B)
  d  <- sqrt(rowSums((pr$A - pr$B)^2))
  
  rms[i] <- sqrt(mean(d^2))
}

per_wing <- data.frame(
  wing_id = wing_ids,
  RMS = rms,
  stringsAsFactors = FALSE
)

write.csv(
  per_wing,
  file.path(out_dir, paste0("digitization_error_per_wing_", species, ".csv")),
  row.names = FALSE
)

# 7.6 Summary -------------------------------------------------------------

summary_df <- data.frame(
  species = species,
  n_wings_replicated = nrow(per_wing),
  mean_RMS   = mean(per_wing$RMS),
  sd_RMS     = sd(per_wing$RMS),
  median_RMS = median(per_wing$RMS),
  q05 = as.numeric(quantile(per_wing$RMS, 0.05)),
  q95 = as.numeric(quantile(per_wing$RMS, 0.95)),
  stringsAsFactors = FALSE
)

write.csv(
  summary_df,
  file.path(out_dir, paste0("digitization_error_summary_", species, ".csv")),
  row.names = FALSE
)

cat("\nDigitization error (Procrustes RMS between replicate digitizations)\n")
cat(
  "Mean RMS:", format(summary_df$mean_RMS, digits = 6),
  " SD:", format(summary_df$sd_RMS, digits = 6),
  " (n =", summary_df$n_wings_replicated, " replicated wings)\n",
  sep = ""
)

cat("Outputs written to:\n  ", out_dir, "\n", sep = "")
