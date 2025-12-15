# 4 - LANDMARK SELECTION --------------------------------------------------
#
# Goal:
#   Empirically determine how many landmarks to retain by
#   progressively removing landmarks with the highest
#   replicate measurement error (ME).
#
# Steps:
#   1) Compute landmark-wise ME (dorsal vs ventral)
#   2) Iteratively remove worst 0–5 landmarks
#   3) For each set:
#        - Mean ME
#        - Replicate-aware FA (bilat.symmetry)
#        - FA / ME ratio
#        - FA stability vs full set (Spearman rho)
#
# Outputs:
#   results/04_landmark_selection_forewings/
#     - landmark_ME_rank.csv
#     - landmark_ME_barplot.png
#     - landmark_set_comparison.csv
#     - FA_ME_ratio_by_set.png

# 4.1 - Libraries ---------------------------------------------------------

if (!requireNamespace("geomorph", quietly = TRUE)) {
  install.packages("geomorph")}

library(geomorph)

# 4.2 - Paths -------------------------------------------------------------

coords_file <- "fw_coords.RDS"
meta_file   <- "fw_metadata.csv"

out_dir <- file.path("results", "4_landmark_selection", "forewings")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 4.3 - Load Data ---------------------------------------------------------

coords_all <- readRDS(coords_file)
meta_all   <- read.csv(meta_file, stringsAsFactors = FALSE)

meta_all$side   <- ifelse(meta_all$side %in% c("L","Left"), "Left", "Right")
meta_all$series <- as.integer(meta_all$series)

# Keep only complete replicate specimens
tab <- with(meta_all, table(specimen_uid, side, series))
uids <- dimnames(tab)[[1]]

keep_uid <- vapply(
  uids,
  function(u) all(tab[u,"Left",] >= 1) && all(tab[u,"Right",] >= 1),
  logical(1)
)

uids <- uids[keep_uid]
keep <- meta_all$specimen_uid %in% uids

coords <- coords_all[,,keep]
meta   <- meta_all[keep,]

cat("Using", length(uids), "specimens with full L/R × dorsal/ventral data\n")

p <- dim(coords)[1]
lm_names <- paste0("LM", seq_len(p))
full_idx <- seq_len(p)

# 4.4 - Landmark Measurement Error ----------------------------------------

compute_landmark_ME <- function(coords, meta) {
  
  gpa <- gpagen(coords, print.progress = FALSE)
  Y <- gpa$coords
  
  key <- paste(meta$specimen_uid, meta$side)
  keys <- unique(key)
  
  dists <- matrix(NA, p, length(keys))
  
  for (i in seq_along(keys)) {
    rows <- which(key == keys[i])
    A <- Y[,, rows[meta$series[rows] == 1][1]]
    B <- Y[,, rows[meta$series[rows] == 2][1]]
    dists[,i] <- sqrt(rowSums((A - B)^2))
  }
  
  data.frame(
    landmark = lm_names,
    ME_mean  = rowMeans(dists),
    stringsAsFactors = FALSE
  )
}

cat("\nComputing landmark measurement error...\n")

me_df <- compute_landmark_ME(coords, meta)
me_df <- me_df[order(me_df$ME_mean, decreasing = TRUE), ]

write.csv(me_df, file.path(out_dir, "landmark_ME_rank.csv"), row.names = FALSE)

png(file.path(out_dir, "landmark_ME_barplot.png"), 1200, 700)
barplot(me_df$ME_mean,
        names.arg = me_df$landmark,
        las = 2,
        main = "Forewing landmark measurement error (series 1 vs 2)",
        ylab = "Mean replicate displacement (Procrustes units)")
dev.off()

# 4.5 - Landmark Set Evaluation -------------------------------------------

run_FA <- function(keep_idx) {
  
  coords_sub <- coords[keep_idx,,,drop=FALSE]
  uids <- sort(unique(meta$specimen_uid))
  psub <- length(keep_idx)
  
  A <- array(NA, c(psub, 2, 4*length(uids)))
  ind <- side <- rep <- character(4*length(uids))
  
  k <- 1
  for (u in uids)
    for (s in c("Left","Right"))
      for (r in c(1,2)) {
        i <- which(meta$specimen_uid==u & meta$side==s & meta$series==r)[1]
        A[,,k] <- coords_sub[,,i]
        ind[k]  <- u
        side[k] <- ifelse(s=="Left","L","R")
        rep[k]  <- r
        k <- k+1
      }
  
  sym <- bilat.symmetry(
    A = A,
    ind = ind,
    side = side,
    replicate = rep,
    object.sym = FALSE,
    RRPP = TRUE,
    print.progress = FALSE
  )
  
  data.frame(
    specimen_uid = names(sym$unsigned.AI),
    FA = as.numeric(sym$unsigned.AI)
  )
}

# Iteratively removing worst landmarks
results <- list()
fa_full <- run_FA(full_idx)

for (k in 0:5) {
  
  removed_lms <- if (k == 0) character(0) else me_df$landmark[1:k]
  keep_idx <- if (k == 0) full_idx else
    setdiff(full_idx, as.integer(sub("LM","", removed_lms)))
  
  fa <- run_FA(keep_idx)
  
  mean_ME <- mean(me_df$ME_mean[keep_idx])
  mean_FA <- mean(fa$FA)
  
  if (k == 0) {
    rho <- 1
  } else {
    merged <- merge(
      fa_full, fa,
      by = "specimen_uid",
      suffixes = c("_full","_red")
    )
    rho <- cor(merged$FA_full, merged$FA_red, method = "spearman")
  }
  
  results[[k+1]] <- data.frame(
    removed_n = k,
    landmarks_used = length(keep_idx),
    removed_landmarks = ifelse(
      length(removed_lms)==0, "none",
      paste(removed_lms, collapse=",")
    ),
    mean_ME = mean_ME,
    mean_FA = mean_FA,
    FA_ME_ratio = mean_FA / mean_ME,
    rho_vs_full = rho,
    stringsAsFactors = FALSE
  )
  
  cat(
    sprintf(
      "Finished set: removed %d landmark(s): %s\n",
      k,
      ifelse(length(removed_lms)==0, "none",
             paste(removed_lms, collapse=", "))
    )
  )
}

res_df <- do.call(rbind, results)
print(res_df)

write.csv(res_df, file.path(out_dir, "landmark_set_comparison.csv"), row.names = FALSE)

# 4.6 - Plot FA / ME Ratio ------------------------------------------------

png(file.path(out_dir, "FA_ME_ratio_by_set.png"), 900, 600)
barplot(
  res_df$FA_ME_ratio,
  names.arg = paste0("-", res_df$removed_n),
  xlab = "Number of worst landmarks removed",
  ylab = "Mean FA / Mean ME",
  main = "Signal-to-noise by landmark set"
)
dev.off()
