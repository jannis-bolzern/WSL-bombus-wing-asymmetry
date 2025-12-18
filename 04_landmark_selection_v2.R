# 04 - LANDMARK SELECTION --------------------------------------------------
#
# Goal:
#   Use replicate photos (series 1 vs 2) to:
#     (1) rank landmarks by replicate measurement error (ME)
#     (2) evaluate a simple fixed removal path (remove worst k landmarks)
#     (3) test a small set of user-defined "awkward" landmarks on top of a baseline
#
# Outputs:
#   results/04_landmark_selection/<wing>/
#     - landmark_ME_rank.csv
#     - landmark_ME_barplot.png
#     - fixed_removal_path.csv
#     - fixed_path_plot.png
#     - awkward_candidates.csv
#     - awkward_candidates_plot.png
#     - recommended_landmark_set.txt

# 4.1 - Libraries ---------------------------------------------------------

if (!requireNamespace("geomorph", quietly = TRUE)) install.packages("geomorph")
library(geomorph)

# 4.2 - Settings ----------------------------------------------------------

wing <- "forewings"  # "forewings" or "hindwings"

coords_file <- if (wing == "forewings") "fw_coords.RDS" else "hw_coords.RDS"
meta_file   <- if (wing == "forewings") "fw_metadata.csv" else "hw_metadata.csv"

out_dir <- file.path("results", "04_landmark_selection", wing)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

max_remove <- 5            # evaluate k = 0..max_remove for fixed removal path
baseline_k <- 2            # baseline = remove top baseline_k landmarks by ME ranking
awkward_lms <- c(6,15,21)  # landmarks to test as awkward candidates (on top of baseline)

# 4.3 - Load data ---------------------------------------------------------

coords_all <- readRDS(coords_file)
meta_all   <- read.csv(meta_file, stringsAsFactors = FALSE)

meta_all$side   <- ifelse(meta_all$side %in% c("L","Left"), "Left", "Right")
meta_all$series <- as.integer(meta_all$series)

# keep only complete L/R × series 1/2 specimens
tab <- with(meta_all, table(specimen_uid, side, series))
uids <- dimnames(tab)[[1]]

complete_uid <- vapply(uids, function(uid) {
  all(
    tab[uid, "Left",  "1"] >= 1,
    tab[uid, "Left",  "2"] >= 1,
    tab[uid, "Right", "1"] >= 1,
    tab[uid, "Right", "2"] >= 1
  )
}, logical(1))

uids_keep <- uids[complete_uid]
sel <- meta_all$specimen_uid %in% uids_keep

coords <- coords_all[,,sel, drop = FALSE]
meta   <- meta_all[sel, , drop = FALSE]

cat("Using", length(uids_keep), "specimens (", nrow(meta), "images) for", wing, "\n")

p <- dim(coords)[1]
full_idx <- seq_len(p)
lm_names <- paste0("LM", full_idx)

max_remove <- min(max_remove, p - 3)
baseline_k <- min(baseline_k, max_remove)

# 4.4 - Helper functions --------------------------------------------------

# Pairwise Procrustes (centroid-size scaling ON by default)
pairwise_procrustes <- function(A, B, scale = TRUE) {
  
  A0 <- sweep(A, 2, colMeans(A), "-")
  B0 <- sweep(B, 2, colMeans(B), "-")
  
  if (scale) {
    csA <- sqrt(sum(A0^2))
    csB <- sqrt(sum(B0^2))
    if (csA > 0) A0 <- A0 / csA
    if (csB > 0) B0 <- B0 / csB
  }
  
  M <- t(B0) %*% A0
  sv <- svd(M)
  R <- sv$u %*% t(sv$v)
  
  if (det(R) < 0) {
    sv$u[, ncol(sv$u)] <- -sv$u[, ncol(sv$u)]
    R <- sv$u %*% t(sv$v)
  }
  
  list(A = A0, B = B0 %*% R)
}

compute_ME <- function(coords_sub, meta_sub) {
  # returns:
  #   ME_vec: per-landmark mean replicate displacement
  #   mean_pair_RMS: overall replicate mismatch (RMS across landmarks, averaged across pairs)
  
  psub <- dim(coords_sub)[1]
  key <- paste(meta_sub$specimen_uid, meta_sub$side)
  keys <- unique(key)
  
  dists <- matrix(NA_real_, nrow = psub, ncol = length(keys))
  
  for (i in seq_along(keys)) {
    rows <- which(key == keys[i])
    i1 <- rows[meta_sub$series[rows] == 1][1]
    i2 <- rows[meta_sub$series[rows] == 2][1]
    
    A <- coords_sub[,,i1]
    B <- coords_sub[,,i2]
    
    pr <- pairwise_procrustes(A, B, scale = TRUE)
    dists[, i] <- sqrt(rowSums((pr$A - pr$B)^2))
  }
  
  list(
    ME_vec = rowMeans(dists, na.rm = TRUE),
    mean_pair_RMS = mean(sqrt(colMeans(dists^2, na.rm = TRUE)), na.rm = TRUE)
  )
}

build_bilat_input <- function(coords_sub, meta_sub) {
  
  uids_sorted <- sort(unique(meta_sub$specimen_uid))
  psub <- dim(coords_sub)[1]
  kdim <- dim(coords_sub)[2]
  
  A <- array(NA_real_, dim = c(psub, kdim, 4 * length(uids_sorted)))
  ind_vec <- character(4 * length(uids_sorted))
  side_vec <- character(4 * length(uids_sorted))
  rep_vec <- integer(4 * length(uids_sorted))
  
  idx <- 1
  for (uid in uids_sorted) {
    for (side in c("Left","Right")) {
      for (rep in c(1,2)) {
        rows <- which(meta_sub$specimen_uid == uid &
                        meta_sub$side == side &
                        meta_sub$series == rep)
        if (length(rows) != 1) stop("Duplicate/missing row for ", uid, " ", side, " rep ", rep)
        A[,,idx] <- coords_sub[,,rows]
        ind_vec[idx] <- uid
        side_vec[idx] <- ifelse(side == "Left", "L", "R")
        rep_vec[idx] <- rep
        idx <- idx + 1
      }
    }
  }
  
  list(
    A = A,
    ind = factor(ind_vec),
    side = factor(side_vec, levels = c("L","R")),
    replicate = factor(rep_vec)
  )
}

run_FA <- function(keep_idx) {
  
  coords_sub <- coords[keep_idx,,, drop = FALSE]
  inp <- build_bilat_input(coords_sub, meta)
  
  sym <- bilat.symmetry(
    A = inp$A,
    ind = inp$ind,
    side = inp$side,
    replicate = inp$replicate,
    object.sym = FALSE,
    RRPP = TRUE,
    iter = 999,
    print.progress = FALSE
  )
  
  data.frame(
    specimen_uid = names(sym$unsigned.AI),
    FA_unsigned = as.numeric(sym$unsigned.AI),
    stringsAsFactors = FALSE
  )
}

rho_between <- function(fa_a, fa_b, suf_a = "_a", suf_b = "_b") {
  m <- merge(fa_a, fa_b, by = "specimen_uid", suffixes = c(suf_a, suf_b))
  cor(m[[paste0("FA_unsigned", suf_a)]], m[[paste0("FA_unsigned", suf_b)]], method = "spearman")
}

make_subsets <- function(v) {
  out <- list(integer(0))
  if (length(v) == 0) return(out)
  for (k in 1:length(v)) out <- c(out, combn(v, k, simplify = FALSE))
  out
}

# 4.5 - Full-set ME rank --------------------------------------------------

cat("\nComputing full-set replicate ME...\n")
ME_full <- compute_ME(coords[full_idx,,, drop = FALSE], meta)

me_rank <- data.frame(
  landmark = lm_names,
  landmark_id = full_idx,
  ME_mean = ME_full$ME_vec,
  stringsAsFactors = FALSE
)
me_rank <- me_rank[order(me_rank$ME_mean, decreasing = TRUE), ]

write.csv(me_rank, file.path(out_dir, "landmark_ME_rank.csv"), row.names = FALSE)

png(file.path(out_dir, "landmark_ME_barplot.png"), 1200, 700)
barplot(me_rank$ME_mean,
        names.arg = me_rank$landmark,
        las = 2,
        main = paste("Landmark replicate ME (pairwise Procrustes) -", wing),
        ylab = "Mean replicate displacement (Procrustes units)")
dev.off()

# 4.6 - Fixed removal path ------------------------------------------------

cat("\nEvaluating fixed removal path...\n")

fa_full <- run_FA(full_idx)

fixed_list <- vector("list", max_remove + 1)

for (k in 0:max_remove) {
  
  removed_ids <- if (k == 0) integer(0) else me_rank$landmark_id[1:k]
  keep_idx <- setdiff(full_idx, removed_ids)
  
  ME_set <- compute_ME(coords[keep_idx,,, drop = FALSE], meta)
  fa_set <- run_FA(keep_idx)
  
  rho <- if (k == 0) 1 else rho_between(fa_full, fa_set, "_full", "_red")
  
  fixed_list[[k + 1]] <- data.frame(
    removed_n = k,
    landmarks_used = length(keep_idx),
    removed_landmarks = if (k == 0) "none" else paste0("LM", removed_ids, collapse = ","),
    mean_pair_RMS = ME_set$mean_pair_RMS,
    rho_vs_full = rho,
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("  removed %d: %s\n", k, fixed_list[[k + 1]]$removed_landmarks))
}

fixed_df <- do.call(rbind, fixed_list)
print(fixed_df)
write.csv(fixed_df, file.path(out_dir, "fixed_removal_path.csv"), row.names = FALSE)

png(file.path(out_dir, "fixed_path_plot.png"), width = 2600, height = 1600, res = 300)
par(mfrow = c(2, 1), mar = c(4, 5, 2.5, 1), oma = c(2, 0, 2, 0))
x <- fixed_df$removed_n

plot(x, fixed_df$rho_vs_full, type = "b", pch = 16, ylim = c(0, 1),
     ylab = expression(paste("Spearman ", rho, " vs full")),
     xlab = "", main = paste("FA stability vs full (fixed removal) -", wing),
     xaxt = "n")
abline(h = 0.95, lty = 2, col = "grey50")
axis(1, at = x)

plot(x, fixed_df$mean_pair_RMS, type = "b", pch = 16,
     ylab = "Mean replicate mismatch",
     xlab = "Number of worst landmarks removed",
     main = "Replicate ME (pairwise Procrustes)")
axis(1, at = x)

dev.off()

# 4.7 - Baseline set ------------------------------------------------------

baseline_removed <- if (baseline_k == 0) integer(0) else me_rank$landmark_id[1:baseline_k]
keep_baseline <- setdiff(full_idx, baseline_removed)

baseline_ME <- compute_ME(coords[keep_baseline,,, drop = FALSE], meta)
baseline_fa <- run_FA(keep_baseline)

cat("\nBaseline: remove ", baseline_k, " landmark(s): ",
    if (baseline_k == 0) "none" else paste0("LM", baseline_removed, collapse = ","),
    "\nBaseline mean_pair_RMS = ", format(baseline_ME$mean_pair_RMS, digits = 6),
    "\n", sep = "")

# 4.8 - Awkward landmark candidates --------------------------------------

cat("\nTesting awkward landmark candidates...\n")

subs <- make_subsets(awkward_lms)

awk_list <- lapply(subs, function(drop_extra) {
  
  keep_idx <- setdiff(keep_baseline, drop_extra)
  
  ME_set <- compute_ME(coords[keep_idx,,, drop = FALSE], meta)
  fa_set <- run_FA(keep_idx)
  
  data.frame(
    drop_extra = if (length(drop_extra) == 0) "none" else paste0("LM", drop_extra, collapse = ","),
    extra_n = length(drop_extra),
    landmarks_used = length(keep_idx),
    mean_pair_RMS = ME_set$mean_pair_RMS,
    rho_vs_baseline = rho_between(baseline_fa, fa_set, "_base", "_new"),
    stringsAsFactors = FALSE
  )
})

awk_df <- do.call(rbind, awk_list)
awk_df <- awk_df[order(awk_df$mean_pair_RMS), ]
print(awk_df)
write.csv(awk_df, file.path(out_dir, "awkward_candidates.csv"), row.names = FALSE)

# choose recommendation
eligible <- awk_df[awk_df$rho_vs_baseline >= 0.95, , drop = FALSE]

eligible <- eligible[eligible$mean_pair_RMS < baseline_ME$mean_pair_RMS, , drop = FALSE]


if (nrow(eligible) == 0) {
  rec_label <- "none"
  rec_removed_extra <- integer(0)
  rec_keep <- keep_baseline
} else {
  eligible <- eligible[order(eligible$mean_pair_RMS, eligible$extra_n), ]
  rec_label <- eligible$drop_extra[1]
  rec_removed_extra <- if (rec_label == "none") integer(0) else as.integer(gsub("LM", "", strsplit(rec_label, ",")[[1]]))
  rec_keep <- setdiff(keep_baseline, rec_removed_extra)
}

rec_removed <- sort(unique(c(baseline_removed, rec_removed_extra)))
rec_keep_ids <- setdiff(full_idx, rec_removed)

writeLines(
  c(
    paste("Wing:", wing),
    paste("Baseline removed:", if (length(baseline_removed) == 0) "none" else paste0("LM", baseline_removed, collapse = ",")),
    paste("Awkward tested:", if (length(awkward_lms) == 0) "none" else paste0("LM", awkward_lms, collapse = ",")),
    paste("rho threshold:", 0.95),
    "",
    paste("RECOMMENDED removed:", if (length(rec_removed) == 0) "none" else paste0("LM", rec_removed, collapse = ",")),
    paste("RECOMMENDED keep:", paste0("LM", rec_keep_ids, collapse = ",")),
    ""
  ),
  con = file.path(out_dir, "recommended_landmark_set.txt")
)

png(file.path(out_dir, "awkward_candidates_plot.png"), width = 2400, height = 1600, res = 300)
par(mar = c(5, 5, 3, 1))

plot(awk_df$mean_pair_RMS, awk_df$rho_vs_baseline,
     pch = 16,
     xlab = "Mean replicate mismatch (pair RMS)",
     ylab = expression(paste("Spearman ", rho, " vs baseline")),
     main = paste("Awkward landmark candidates (relative to baseline) -", wing),
     ylim = c(min(0.9, min(awk_df$rho_vs_baseline, na.rm = TRUE)), 1.0))
abline(h = 0.95, lty = 2, col = "grey50")

text(awk_df$mean_pair_RMS, awk_df$rho_vs_baseline, labels = awk_df$drop_extra, pos = 3, cex = 0.7)

# baseline point
points(baseline_ME$mean_pair_RMS, 1, pch = 17, cex = 1.3)

# recommended point
rec_ME <- compute_ME(coords[rec_keep,,, drop = FALSE], meta)$mean_pair_RMS
rec_rho <- rho_between(baseline_fa, run_FA(rec_keep), "_base", "_new")
points(rec_ME, rec_rho, pch = 8, cex = 1.6)

legend("bottomleft", legend = c("Candidates", "Baseline", "Recommended"),
       pch = c(16, 17, 8), bty = "n", cex = 0.9)

dev.off()

cat("\nDone. Outputs written to:\n  ", out_dir, "\n", sep = "")
