# 05 - FLUCTUATING ASYMMETRY EXTRACTION -----------------------------------
#
# Purpose:
#   - Extract fluctuating asymmetry (FA) from forewing shape
#   - Single digitization per wing (series = 1)
#   - Replicates used only for QC, not FA modeling
#   - Run separately for Bombus lapidarius and B. pascuorum
#
# Outputs (per species):
#   results/05_fa_analysis/
#     - fa_table_<SPECIES>.csv
#     - symmetry_object_<SPECIES>.RDS

# 5.1 - Libraries ---------------------------------------------------------

if (!requireNamespace("geomorph", quietly = TRUE)) {
  install.packages("geomorph")
}
library(geomorph)

# 5.2 - Settings ----------------------------------------------------------

data_dir <- "data"
coords_prefix <- "fw_coords_"
meta_prefix   <- "fw_metadata_"

out_dir <- file.path("results", "05_extract_fa")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 5.3 - Detect species automatically -------------------------------------

coord_files <- list.files(
  data_dir,
  pattern = paste0("^", coords_prefix, ".*\\.RDS$"),
  full.names = FALSE
)

if (length(coord_files) == 0) {
  stop("No fw_coords_<species>.RDS files found in ", data_dir)
}

species <- gsub(
  paste0("^", coords_prefix, "|\\.RDS$"),
  "",
  coord_files
)

cat("Species detected for FA analysis:\n")
for (sp in species) cat("  - ", sp, "\n", sep = "")

# 5.4 - Loop over species -------------------------------------------------
for (sp in species) {
  
  coords_file <- file.path(data_dir, paste0(coords_prefix, sp, ".RDS"))
  meta_file   <- file.path(data_dir, paste0(meta_prefix, sp, ".csv"))
  
  if (!file.exists(coords_file) || !file.exists(meta_file)) {
    message("Skipping species ", sp, " (missing data files).")
    next
  }
  
  cat("\nProcessing species:", sp, "\n")
  
  coords_all <- readRDS(coords_file)
  meta_all   <- read.csv(meta_file, stringsAsFactors = FALSE)
  
  if (nrow(meta_all) != dim(coords_all)[3]) {
    stop("Metadata / coordinate mismatch for species: ", sp)
  }
  
  # Normalize side labels
  meta_all$side <- ifelse(meta_all$side %in% c("L", "Left"), "Left", "Right")
  
  # Filter to primary digitization only
  keep <- meta_all$series == 1
  
  cat(
    "  Dropping", sum(!keep),
    "replicate digitizations\n"
  )
  
  meta_all   <- meta_all[keep, , drop = FALSE]
  coords_all <- coords_all[,,keep, drop = FALSE]
  
  # Check for Left / Right availability
  if (!all(c("Left", "Right") %in% meta_all$side)) {
    stop("Missing Left or Right wings for species: ", sp)
  }
  
  tab <- with(meta_all, table(specimen_uid, side))
  
  uids_keep <- rownames(tab)[
    tab[, "Left"] == 1 &
      tab[, "Right"] == 1
  ]
  
  if (length(uids_keep) == 0) {
    stop("No complete Left/Right specimens for species: ", sp)
  }
  
  sel <- meta_all$specimen_uid %in% uids_keep
  coords <- coords_all[,,sel, drop = FALSE]
  meta   <- meta_all[sel, , drop = FALSE]
  
  cat("  Specimens used:", length(uids_keep), "\n")
  cat("  Wings used:", nrow(meta), "\n")
  
  # Build array for bilat.symmetry
  uids_sorted <- sort(unique(meta$specimen_uid))
  p <- dim(coords)[1]
  k <- dim(coords)[2]
  
  A <- array(NA_real_, dim = c(p, k, 2 * length(uids_sorted)))
  ind_vec  <- character(2 * length(uids_sorted))
  side_vec <- character(2 * length(uids_sorted))
  
  idx <- 1
  for (uid in uids_sorted) {
    for (side in c("Left", "Right")) {
      
      rows <- which(meta$specimen_uid == uid & meta$side == side)
      
      if (length(rows) != 1) {
        stop("Unexpected number of wings for ", uid, " (", side, ")")
      }
      
      A[,,idx] <- coords[,,rows]
      ind_vec[idx]  <- uid
      side_vec[idx] <- ifelse(side == "Left", "L", "R")
      idx <- idx + 1
    }
  }
  
  # Run FA
  sym <- bilat.symmetry(
    A = A,
    ind = factor(ind_vec),
    side = factor(side_vec, levels = c("L","R")),
    object.sym = FALSE,
    RRPP = TRUE,
    iter = 999,
    print.progress = FALSE
  )
  
  # Extract FA values
  fa <- data.frame(
    specimen_uid = names(sym$unsigned.AI),
    FA_unsigned  = as.numeric(sym$unsigned.AI),
    stringsAsFactors = FALSE
  )
  
  spec_meta <- unique(meta[, c("specimen_uid", "urbanity", "species", "site")])
  fa <- merge(fa, spec_meta, by = "specimen_uid", all.x = TRUE)
  
  # Write outputs
  write.csv(
    fa,
    file.path(out_dir, paste0("fa_table_", sp, ".csv")),
    row.names = FALSE
  )
  
  saveRDS(
    sym,
    file.path(out_dir, paste0("symmetry_object_", sp, ".RDS"))
  )
  
  cat("  FA extraction complete for", sp, "\n")
}

cat("\nOutputs written to:\n  ", out_dir, "\n", sep = "")
