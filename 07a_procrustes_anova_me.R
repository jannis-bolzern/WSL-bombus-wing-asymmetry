# 07a - PROCRUSTES ANOVA: MEASUREMENT ERROR VALIDATION --------------------
#
# Runs TWO clean ANOVAs per species:
#   A) Digitization ME: photo_rep==1, digit_rep 1 vs 2
#   B) Imaging ME: digit_rep==1, photo_rep 1 vs 2
#
# Outputs (per species):
#   results/07a_procrustes_anova/
#     - procrustes_anova_digit_ME_summary_<SPECIES>.txt
#     - procrustes_anova_digit_ME_table_<SPECIES>.csv
#     - procrustes_anova_image_ME_summary_<SPECIES>.txt
#     - procrustes_anova_image_ME_table_<SPECIES>.csv
#     - procrustes_anova_ME_table_<SPECIES>_combined.csv   (NEW)

# 7a.1 - Libraries --------------------------------------------------------

if (!requireNamespace("geomorph", quietly = TRUE)) install.packages("geomorph")
library(geomorph)

# 7a.2 - Settings ---------------------------------------------------------

data_dir <- "data"
coords_prefix <- "fw_coords_"
meta_prefix   <- "fw_metadata_"

out_dir <- file.path("results", "07a_procrustes_anova")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

coord_files <- list.files(
  data_dir,
  pattern = paste0("^", coords_prefix, ".*\\.RDS$"),
  full.names = FALSE
)
if (length(coord_files) == 0) stop("No fw_coords_<species>.RDS found in ", data_dir)

species <- tools::file_path_sans_ext(sub(paste0("^", coords_prefix), "", coord_files))

cat("Species detected for Procrustes ANOVA:\n")
for (sp in species) cat("  - ", sp, "\n", sep = "")

# 7a.3 - Helpers ----------------------------------------------------------

# Robustly extract the ANOVA table across geomorph/RRPP versions
extract_aov_table <- function(fit) {
  sm <- tryCatch(summary(fit), error = function(e) NULL)
  
  candidates <- list(
    if (!is.null(sm)) sm$aov.table else NULL,
    if (!is.null(sm)) sm$table else NULL,
    fit$aov.table,
    fit$anova.table,
    fit$ANOVA
  )
  
  tab <- NULL
  for (x in candidates) {
    if (!is.null(x)) { tab <- x; break }
  }
  if (is.null(tab)) tab <- tryCatch(anova(fit), error = function(e) NULL)
  if (is.null(tab)) stop("Could not extract ANOVA table from procD.lm object.")
  
  df <- as.data.frame(tab)
  df$term_raw <- rownames(df)
  rownames(df) <- NULL
  df
}

write_summary_txt <- function(fit, file, header_lines = character(0)) {
  sm <- summary(fit)
  sink(file)
  if (length(header_lines) > 0) cat(paste0(header_lines, collapse = "\n"), "\n\n")
  print(sm)
  sink()
}

# publication-ready formatting
sig_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  if (p < 0.1)   return(".")
  ""
}

clean_term <- function(term_raw, test_type) {
  # term names depend on how procD.lm labels your formula terms
  # We standardize the common ones:
  if (term_raw %in% c("ind")) return("Individual")
  if (term_raw %in% c("sideF")) return("Side (DA)")
  if (term_raw %in% c("ind:sideF")) return("Individual × Side (FA)")
  if (term_raw %in% c("digrep")) return("Digitization replicate")
  if (term_raw %in% c("photorep")) return("Photo replicate")
  if (grepl("Residual", term_raw, ignore.case = TRUE)) return("Residuals")
  # fallback:
  term_raw
}

# find p-value column robustly
get_p_col <- function(df) {
  p_candidates <- names(df)[grepl("^Pr\\(", names(df)) | grepl("p", names(df), ignore.case = TRUE)]
  # prefer "Pr(>F)" if it exists
  if ("Pr(>F)" %in% names(df)) return("Pr(>F)")
  if (length(p_candidates) > 0) return(p_candidates[1])
  NULL
}

# build publication-ready table
format_pub_table <- function(df_raw, sp, test_label, n_ind, n_wings, n_obs) {
  
  # Map columns robustly
  df <- df_raw
  
  # df
  if (!"Df" %in% names(df)) {
    if ("df" %in% names(df)) df$Df <- df$df else df$Df <- NA_integer_
  }
  # SS, MS
  if (!"SS" %in% names(df)) df$SS <- NA_real_
  if (!"MS" %in% names(df)) df$MS <- NA_real_
  # R2 (geomorph uses Rsq)
  if (!"Rsq" %in% names(df)) {
    if ("R2" %in% names(df)) df$Rsq <- df$R2 else df$Rsq <- NA_real_
  }
  # F, Z
  if (!"F" %in% names(df)) df$F <- NA_real_
  if (!"Z" %in% names(df)) df$Z <- NA_real_
  
  pcol <- get_p_col(df)
  df$p <- if (!is.null(pcol)) df[[pcol]] else NA_real_
  
  # clean terms
  df$Source <- vapply(df$term_raw, clean_term, character(1), test_type = test_label)
  
  # order (typical Procrustes ANOVA order)
  ord_levels <- c("Individual", "Side (DA)", "Individual × Side (FA)",
                  "Digitization replicate", "Photo replicate", "Residuals")
  df$Source <- factor(df$Source, levels = ord_levels, ordered = TRUE)
  
  # build final
  out <- data.frame(
    Species = sp,
    Test = test_label,
    N_ind = n_ind,
    N_wings = n_wings,
    N_obs = n_obs,
    
    Source = as.character(df$Source),
    Df = as.integer(df$Df),
    
    SS = as.numeric(df$SS),
    MS = as.numeric(df$MS),
    R2 = as.numeric(df$Rsq),
    R2_pct = 100 * as.numeric(df$Rsq),
    
    F = as.numeric(df$F),
    Z = as.numeric(df$Z),
    p = as.numeric(df$p),
    
    sig = vapply(as.numeric(df$p), sig_stars, character(1)),
    stringsAsFactors = FALSE
  )
  
  # drop rows with NA source (shouldn't happen, but safe)
  out <- out[!is.na(out$Source) & out$Source != "", , drop = FALSE]
  
  # sort by Source factor order (Residuals last)
  out$Source <- factor(out$Source, levels = ord_levels, ordered = TRUE)
  out <- out[order(out$Source), , drop = FALSE]
  out$Source <- as.character(out$Source)
  
  # Nice rounding for publication tables (still numeric in CSV)
  out$SS <- signif(out$SS, 6)
  out$MS <- signif(out$MS, 6)
  out$R2 <- round(out$R2, 6)
  out$R2_pct <- round(out$R2_pct, 4)
  out$F <- round(out$F, 3)
  out$Z <- round(out$Z, 3)
  out$p <- signif(out$p, 4)
  
  out
}

# 7a.4 - Loop over species -------------------------------------------------

for (sp in species) {
  
  coords_file <- file.path(data_dir, paste0(coords_prefix, sp, ".RDS"))
  meta_file   <- file.path(data_dir, paste0(meta_prefix, sp, ".csv"))
  
  if (!file.exists(coords_file) || !file.exists(meta_file)) {
    message("Skipping species ", sp, " (missing files).")
    next
  }
  
  cat("\nProcessing species:", sp, "\n")
  
  coords_all <- readRDS(coords_file)
  meta_all   <- read.csv(meta_file, stringsAsFactors = FALSE)
  
  if (nrow(meta_all) != dim(coords_all)[3]) {
    stop("Metadata / coordinate mismatch for species: ", sp)
  }
  
  need_cols <- c("specimen_uid", "side", "photo_rep", "digit_rep")
  if (!all(need_cols %in% names(meta_all))) {
    stop("Missing required columns in metadata for ", sp, ": ",
         paste(setdiff(need_cols, names(meta_all)), collapse = ", "))
  }
  
  meta_all$side <- ifelse(meta_all$side %in% c("L","Left"), "L", "R")
  meta_all$photo_rep <- as.integer(meta_all$photo_rep)
  meta_all$digit_rep <- as.integer(meta_all$digit_rep)
  
  rep_uids <- unique(meta_all$specimen_uid[meta_all$photo_rep == 2])
  if (length(rep_uids) == 0) {
    message("  No imaging replicate individuals found for ", sp, " (photo_rep==2).")
    next
  }
  
  combined_tables <- list()
  
  # ----------------------------------------------------------------------
  # A) DIGITIZATION ME
  # ----------------------------------------------------------------------
  
  idx_dig <- which(meta_all$specimen_uid %in% rep_uids &
                     meta_all$photo_rep == 1 &
                     meta_all$digit_rep %in% c(1, 2))
  
  meta_dig <- meta_all[idx_dig, , drop = FALSE]
  coords_dig <- coords_all[,, idx_dig, drop = FALSE]
  meta_dig$wing_id <- paste(meta_dig$specimen_uid, meta_dig$side, sep = "__")
  
  tab_dig <- with(meta_dig, table(wing_id, digit_rep))
  if (all(c("1","2") %in% colnames(tab_dig))) {
    
    keep_wings <- rownames(tab_dig)[tab_dig[, "1"] == 1 & tab_dig[, "2"] == 1]
    keep_idx <- meta_dig$wing_id %in% keep_wings
    
    meta_dig <- meta_dig[keep_idx, , drop = FALSE]
    coords_dig <- coords_dig[,, keep_idx, drop = FALSE]
    
    if (nrow(meta_dig) > 0) {
      
      gpa <- gpagen(coords_dig, print.progress = FALSE)
      
      meta_dig$ind   <- factor(meta_dig$specimen_uid)
      meta_dig$sideF <- factor(meta_dig$side, levels = c("L","R"))
      meta_dig$digrep <- factor(meta_dig$digit_rep)
      
      fit_dig <- procD.lm(gpa$coords ~ ind * sideF + digrep, data = meta_dig, iter = 999)
      
      sum_file <- file.path(out_dir, paste0("procrustes_anova_digit_ME_summary_", sp, ".txt"))
      write_summary_txt(
        fit_dig, sum_file,
        header_lines = c(paste("Species:", sp),
                         "Test: Digitization ME (photo_rep==1, digit_rep 1 vs 2)")
      )
      
      raw_tab <- extract_aov_table(fit_dig)
      
      n_ind <- length(unique(meta_dig$specimen_uid))
      n_wings <- length(unique(meta_dig$wing_id))
      n_obs <- nrow(meta_dig)
      
      pub_tab <- format_pub_table(raw_tab, sp, "Digitization ME", n_ind, n_wings, n_obs)
      
      out_csv <- file.path(out_dir, paste0("procrustes_anova_digit_ME_table_", sp, ".csv"))
      write.csv(pub_tab, out_csv, row.names = FALSE)
      
      combined_tables[[length(combined_tables) + 1]] <- pub_tab
      
      cat("  Saved digitization ME table for ", sp, "\n", sep = "")
    }
  }
  
  # ----------------------------------------------------------------------
  # B) IMAGING ME
  # ----------------------------------------------------------------------
  
  idx_img <- which(meta_all$specimen_uid %in% rep_uids &
                     meta_all$digit_rep == 1 &
                     meta_all$photo_rep %in% c(1, 2))
  
  meta_img <- meta_all[idx_img, , drop = FALSE]
  coords_img <- coords_all[,, idx_img, drop = FALSE]
  meta_img$wing_id <- paste(meta_img$specimen_uid, meta_img$side, sep = "__")
  
  tab_img <- with(meta_img, table(wing_id, photo_rep))
  if (all(c("1","2") %in% colnames(tab_img))) {
    
    keep_wings <- rownames(tab_img)[tab_img[, "1"] == 1 & tab_img[, "2"] == 1]
    keep_idx <- meta_img$wing_id %in% keep_wings
    
    meta_img <- meta_img[keep_idx, , drop = FALSE]
    coords_img <- coords_img[,, keep_idx, drop = FALSE]
    
    if (nrow(meta_img) > 0) {
      
      gpa <- gpagen(coords_img, print.progress = FALSE)
      
      meta_img$ind   <- factor(meta_img$specimen_uid)
      meta_img$sideF <- factor(meta_img$side, levels = c("L","R"))
      meta_img$photorep <- factor(meta_img$photo_rep)
      
      fit_img <- procD.lm(gpa$coords ~ ind * sideF + photorep, data = meta_img, iter = 999)
      
      sum_file <- file.path(out_dir, paste0("procrustes_anova_image_ME_summary_", sp, ".txt"))
      write_summary_txt(
        fit_img, sum_file,
        header_lines = c(paste("Species:", sp),
                         "Test: Imaging ME (digit_rep==1, photo_rep 1 vs 2)")
      )
      
      raw_tab <- extract_aov_table(fit_img)
      
      n_ind <- length(unique(meta_img$specimen_uid))
      n_wings <- length(unique(meta_img$wing_id))
      n_obs <- nrow(meta_img)
      
      pub_tab <- format_pub_table(raw_tab, sp, "Imaging ME", n_ind, n_wings, n_obs)
      
      out_csv <- file.path(out_dir, paste0("procrustes_anova_image_ME_table_", sp, ".csv"))
      write.csv(pub_tab, out_csv, row.names = FALSE)
      
      combined_tables[[length(combined_tables) + 1]] <- pub_tab
      
      cat("  Saved imaging ME table for ", sp, "\n", sep = "")
    }
  }
  
  # Combined file
  if (length(combined_tables) > 0) {
    comb <- do.call(rbind, combined_tables)
    comb_file <- file.path(out_dir, paste0("procrustes_anova_ME_table_", sp, "_combined.csv"))
    write.csv(comb, comb_file, row.names = FALSE)
    cat("  Saved combined ME table: ", comb_file, "\n", sep = "")
  }
}

cat("\nDone. Outputs in: ", out_dir, "\n", sep = "")
