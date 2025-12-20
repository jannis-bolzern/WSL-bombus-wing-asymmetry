# 06 - LANDMARK SET SENSITIVITY ANALYSIS ----------------------------------
#
# Purpose:
#   - Test robustness of FA estimates to landmark set size
#   - Compare full (21 landmarks) vs reduced (15 landmarks)
#   - Uses forewing FA (single digitization per wing)
#
# Outputs:
#   results/06_landmark_sensitivity/
#     - fa_comparison_table_<SPECIES>.csv
#     - fa_correlation_<SPECIES>.txt
#     - fa_scatterplot_<SPECIES>.png
# ------------------------------------------------------------------------

# 6.1 - Libraries ---------------------------------------------------------

if (!requireNamespace("geomorph", quietly = TRUE)) {
  install.packages("geomorph")
}
library(geomorph)

# 6.2 - Settings ----------------------------------------------------------

species <- "lapidarius"   # <- CHANGE HERE if needed

data_dir <- "data"
coords_file <- file.path(data_dir, paste0("fw_coords_", species, ".RDS"))
meta_file   <- file.path(data_dir, paste0("fw_metadata_", species, ".csv"))

out_dir <- file.path("results", "06_landmark_sensitivity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# landmark sets
lm_full    <- 1:21
lm_reduced <- 1:15

# 6.3 - Load data ---------------------------------------------------------

if (!file.exists(coords_file) || !file.exists(meta_file)) {
  stop("Missing input files for species: ", species)
}

coords_all <- readRDS(coords_file)
meta_all   <- read.csv(meta_file, stringsAsFactors = FALSE)

if (nrow(meta_all) != dim(coords_all)[3]) {
  stop("Metadata / coordinate mismatch for species: ", species)
}

# normalize side labels
meta_all$side <- ifelse(meta_all$side %in% c("L","Left"), "Left", "Right")

# keep only first digitization
keep <- meta_all$series == 1
coords_all <- coords_all[,,keep, drop = FALSE]
meta_all   <- meta_all[keep, , drop = FALSE]

# 6.4 - Helper function: FA for a landmark subset -------------------------

run_fa_subset <- function(coords, meta, keep_lm) {
  
  coords_sub <- coords[keep_lm, , , drop = FALSE]
  
  # keep only complete Left/Right specimens
  tab <- with(meta, table(specimen_uid, side))
  
  if (!all(c("Left","Right") %in% colnames(tab))) {
    stop("Missing Left or Right wings after filtering.")
  }
  
  uids_keep <- rownames(tab)[
    tab[, "Left"] == 1 &
      tab[, "Right"] == 1
  ]
  
  sel <- meta$specimen_uid %in% uids_keep
  coords_sub <- coords_sub[,,sel, drop = FALSE]
  meta_sub   <- meta[sel, , drop = FALSE]
  
  uids_sorted <- sort(unique(meta_sub$specimen_uid))
  p <- dim(coords_sub)[1]
  
  A <- array(NA_real_, dim = c(p, 2, 2 * length(uids_sorted)))
  ind_vec  <- character(2 * length(uids_sorted))
  side_vec <- character(2 * length(uids_sorted))
  
  idx <- 1
  for (uid in uids_sorted) {
    for (side in c("Left","Right")) {
      
      rows <- which(
        meta_sub$specimen_uid == uid &
          meta_sub$side == side
      )
      
      if (length(rows) != 1) {
        stop("Unexpected number of wings for ", uid, " (", side, ")")
      }
      
      A[,,idx] <- coords_sub[,,rows]
      ind_vec[idx]  <- uid
      side_vec[idx] <- ifelse(side == "Left", "L", "R")
      idx <- idx + 1
    }
  }
  
  sym <- bilat.symmetry(
    A = A,
    ind = factor(ind_vec),
    side = factor(side_vec, levels = c("L","R")),
    object.sym = FALSE,
    RRPP = TRUE,
    iter = 999,
    print.progress = FALSE
  )
  
  data.frame(
    specimen_uid = names(sym$unsigned.AI),
    FA_unsigned  = as.numeric(sym$unsigned.AI),
    stringsAsFactors = FALSE
  )
}

# 6.5 - Run FA for both landmark sets -------------------------------------

cat("Running FA with 21 landmarks...\n")
fa_21 <- run_fa_subset(coords_all, meta_all, lm_full)

cat("Running FA with 15 landmarks...\n")
fa_15 <- run_fa_subset(coords_all, meta_all, lm_reduced)

fa_cmp <- merge(
  fa_21, fa_15,
  by = "specimen_uid",
  suffixes = c("_21lm", "_15lm")
)

# 6.6 - Correlation -------------------------------------------------------

rho <- cor(
  fa_cmp$FA_unsigned_21lm,
  fa_cmp$FA_unsigned_15lm,
  method = "spearman"
)

writeLines(
  paste("Spearman correlation (21 vs 15 landmarks):", round(rho, 4)),
  con = file.path(out_dir, paste0("fa_correlation_", species, ".txt"))
)

# 6.7 - Scatterplot -------------------------------------------------------

png(
  file.path(out_dir, paste0("fa_scatterplot_", species, ".png")),
  width = 1600, height = 1600, res = 300
)

plot(
  fa_cmp$FA_unsigned_21lm,
  fa_cmp$FA_unsigned_15lm,
  pch = 16,
  xlab = "FA (21 landmarks)",
  ylab = "FA (15 landmarks)",
  main = paste("FA robustness to landmark reduction\n", species)
)
abline(0, 1, col = "red", lwd = 2)

dev.off()

# 6.8 - Save comparison table --------------------------------------------

write.csv(
  fa_cmp,
  file.path(out_dir, paste0("fa_comparison_table_", species, ".csv")),
  row.names = FALSE
)

# 6.9 - Summary -----------------------------------------------------------

cat("\nSpecies:", species, "\n")
cat("Specimens used:", nrow(fa_cmp), "\n")
cat("Spearman rho:", round(rho, 4), "\n")
cat("Outputs written to:\n  ", out_dir, "\n", sep = "")
