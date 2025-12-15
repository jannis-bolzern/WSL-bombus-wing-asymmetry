# 5 - FA EXTRACTION -------------------------------------------------------
#
# Purpose:
#   - Extract replicate-aware fluctuating asymmetry (FA)
#   - Using a fixed landmark set (chosen in Script 04)
#   - No hypothesis testing here
#
# Output:
#   results/05_FA/FA_per_specimen.csv
#   results/05_FA/symmetry_object.RDS

# 5.1 - Libraries ---------------------------------------------------------

library(geomorph)

# 5.2 - Paths -------------------------------------------------------------

coords_file <- "fw_coords.RDS"
meta_file   <- "fw_metadata.csv"
out_dir <- file.path("results", "5_fa_analysis")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 5.3 - Load Data ---------------------------------------------------------

coords_all <- readRDS(coords_file)
meta_all   <- read.csv(meta_file, stringsAsFactors = FALSE)

meta_all$side   <- ifelse(meta_all$side %in% c("L","Left"), "Left", "Right")
meta_all$series <- as.integer(meta_all$series)

# Keep only complete specimens
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

cat("Using", length(uids_keep), "specimens (",
    nrow(meta), "images ) for FA analysis\n")

# 5.4 - Build Array For bilat.symmetry ------------------------------------

# Order: for each specimen → L1 L2 R1 R2
uids_sorted <- sort(unique(meta$specimen_uid))
p <- dim(coords)[1]
k <- dim(coords)[2]

coords_3d <- array(NA_real_, dim = c(p, k, 4 * length(uids_sorted)))
ind_vec   <- character(4 * length(uids_sorted))
side_vec  <- character(4 * length(uids_sorted))
rep_vec   <- integer(4 * length(uids_sorted))

idx <- 1
for (uid in uids_sorted) {
  for (side in c("Left","Right")) {
    for (rep in c(1,2)) {
      row <- which(meta$specimen_uid == uid &
                     meta$side == side &
                     meta$series == rep)[1]
      
      coords_3d[,,idx] <- coords[,,row]
      ind_vec[idx]  <- uid
      side_vec[idx] <- ifelse(side == "Left", "L", "R")
      rep_vec[idx]  <- rep
      idx <- idx + 1
    }
  }
}

# 5.5 - Run FA (replicate-aware) ------------------------------------------

sym <- bilat.symmetry(
  A = coords_3d,
  ind = factor(ind_vec),
  side = factor(side_vec, levels = c("L","R")),
  replicate = factor(rep_vec),
  object.sym = FALSE,
  RRPP = TRUE,
  iter = 999,
  print.progress = TRUE
)

# 5.6 - Extract FA Table --------------------------------------------------

fa <- data.frame(
  specimen_uid = names(sym$unsigned.AI),
  FA_unsigned  = as.numeric(sym$unsigned.AI),
  stringsAsFactors = FALSE
)

spec_meta <- unique(meta[, c("specimen_uid", "urbanity", "species", "site")])
fa <- merge(fa, spec_meta, by = "specimen_uid", all.x = TRUE)

print(fa)

# 5.7 - Save Outputs ------------------------------------------------------

write.csv(fa,
          file.path(out_dir, "05_forewing_fa_table.csv"),
          row.names = FALSE)

saveRDS(sym,
        file.path(out_dir, "05_forewing_symmetry_object.RDS"))

cat("FA extraction complete.\n")
cat("Saved:\n")
cat(" - results/05_forewing_fa_table.csv\n")
cat(" - results/05_forewing_symmetry_object.RDS\n")
