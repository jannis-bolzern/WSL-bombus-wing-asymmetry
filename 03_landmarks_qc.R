# 03 - QC: LANDMARK COUNT AND ORDER
#
# Philosophy:
#   - FAIL when landmark COUNT is wrong
#   - FAIL when landmark LABEL ORDER is likely wrong
#   - IGNORE pure shape outliers with no swap improvement

# 3.1 - Libraries ---------------------------------------------------------

if (!requireNamespace("StereoMorph", quietly = TRUE))
  install.packages("StereoMorph")
if (!requireNamespace("geomorph", quietly = TRUE))
  install.packages("geomorph")

library(StereoMorph)
library(geomorph)

# 3.2 - Landmark expectations ---------------------------------------------

landmark_spec <- list(
  forewings = 21,
  hindwings = 6
)

base_data_dir <- "data"

# 3.3 - QC 1: Landmark count ----------------------------------------------

qc_check_landmark_count <- function(folder, expected_n) {
  
  files <- list.files(folder, "\\.txt$", full.names = TRUE)
  if (!length(files)) return(invisible(NULL))
  
  bad <- character(0)
  
  for (f in files) {
    shp <- tryCatch(readShapes(f), error = function(e) NULL)
    lm  <- if (!is.null(shp)) shp$landmarks.pixel else NULL
    
    if (is.null(lm) || nrow(lm) != expected_n || ncol(lm) != 2)
      bad <- c(bad, basename(f))
  }
  
  if (length(bad)) {
    cat("\nQC FAILED (landmark count) in", folder, "\n")
    for (x in bad) cat("  -", x, "\n")
    stop("Fix landmark count errors.")
  }
  
  cat("QC PASSED (count):", folder, "\n")
}

# 3.4 - QC 2: Landmark order ----------------------------------------------

qc_check_landmark_order <- function(
    folder,
    expected_n,
    mad_cutoff  = 4.5,
    min_improve = 0.12
) {
  
  files <- list.files(folder, "\\.txt$", full.names = TRUE)
  mats <- lapply(files, function(f) {
    shp <- tryCatch(readShapes(f), error = function(e) NULL)
    if (is.null(shp)) return(NULL)
    as.matrix(shp$landmarks.pixel)
  })
  
  ok <- vapply(
    mats,
    function(m) !is.null(m) && nrow(m) == expected_n && ncol(m) == 2,
    logical(1)
  )
  
  mats  <- mats[ok]
  files <- files[ok]
  
  if (length(mats) < 5) return(invisible(TRUE))
  
  # GPA
  arr <- array(NA_real_, c(expected_n, 2, length(mats)))
  for (i in seq_along(mats)) arr[, , i] <- mats[[i]]
  
  gpa <- gpagen(arr, print.progress = FALSE)
  consensus <- gpa$consensus
  coords <- gpa$coords
  
  # Distances
  d <- vapply(
    seq_len(dim(coords)[3]),
    function(i) sqrt(sum((coords[, , i] - consensus)^2)),
    numeric(1)
  )
  
  mad_z <- (d - median(d)) / mad(d, constant = 1)
  out   <- which(mad_z > mad_cutoff)
  
  if (!length(out)) {
    cat("QC PASSED (order):", folder, "\n")
    return(invisible(TRUE))
  }
  
  likely_swaps <- list()
  
  for (ii in out) {
    
    orig   <- coords[, , ii]
    orig_d <- d[ii]
    
    # Landmark contributions
    lm_disp   <- rowSums((orig - consensus)^2)
    worst_idx <- order(lm_disp, decreasing = TRUE)[1:3]
    
    # Try all swaps
    best_d    <- orig_d
    best_pair <- c(NA_integer_, NA_integer_)
    
    for (p in combn(expected_n, 2, simplify = FALSE)) {
      tmp <- orig
      tmp[p[1], ] <- orig[p[2], ]
      tmp[p[2], ] <- orig[p[1], ]
      td <- sqrt(sum((tmp - consensus)^2))
      
      if (td < best_d) {
        best_d <- td
        best_pair <- p
      }
    }
    
    improve <- (orig_d - best_d) / orig_d
    
    # Require improvement AND localization
    swap_supported <- all(best_pair %in% worst_idx)
    
    if (is.finite(improve) &&
        improve >= min_improve &&
        swap_supported) {
      
      likely_swaps[[length(likely_swaps) + 1]] <- data.frame(
        file = basename(files[ii]),
        dist = orig_d,
        mad_z = mad_z[ii],
        worst_landmarks =
          paste0("LM", worst_idx, collapse = ","),
        suggested_swap =
          paste0("LM", best_pair[1], "<->LM", best_pair[2]),
        improve_pct = 100 * improve,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Decision
  if (length(likely_swaps) == 0) {
    cat("QC PASSED (order):", folder, "\n")
    return(invisible(TRUE))
  }
  
  cat("\nORDER-QC FAILED (LIKELY LANDMARK SWAPS) in", folder, "\n")
  
  df <- do.call(rbind, likely_swaps)
  df <- df[order(-df$improve_pct), ]
  
  print(within(df, {
    dist <- round(dist, 4)
    mad_z <- round(mad_z, 2)
    improve_pct <- round(improve_pct, 0)
  }), row.names = FALSE)
  
  stop("Order-QC failed: likely landmark label swaps detected.")
}

# -------------------------------------------------------------------------
# Run QC
# -------------------------------------------------------------------------

cat("\n--- Running landmark QC ---\n\n")

for (wing in names(landmark_spec)) {
  
  wing_dir <- file.path(base_data_dir, paste0("landmark_data_", wing))
  if (!dir.exists(wing_dir)) next
  
  for (folder in list.dirs(wing_dir, recursive = FALSE)) {
    
    cat("\nChecking", basename(folder), "-", wing, "\n")
    
    qc_check_landmark_count(folder, landmark_spec[[wing]])
    qc_check_landmark_order(folder, landmark_spec[[wing]])
  }
}
