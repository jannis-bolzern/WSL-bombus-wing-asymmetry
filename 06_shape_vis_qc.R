# ============================================================
# 05B — SHAPE VISUALIZATION & QC (FOREWINGS)
# ============================================================
# Produces:
#   - Consensus landmark configuration (after GPA)
#   - PCA of symmetric shape
#   - PCA of asymmetric shape (FA component)
#   - PCA-based QC table (potential outliers)
#
# Outputs saved to:
#   results/05B_shape_plots/
# ============================================================

library(geomorph)

# ---------------------------
# Paths
# ---------------------------
coords_file <- "fw_coords.RDS"
meta_file   <- "fw_metadata.csv"
sym_file    <- "results/05_FA/05_forewing_symmetry_object.RDS"

out_dir  <- "results/05B_shape_plots"
qc_dir   <- file.path(out_dir, "QC")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Load data
# ---------------------------
coords <- readRDS(coords_file)        # p x 2 x n
meta   <- read.csv(meta_file, stringsAsFactors = FALSE)
sym    <- readRDS(sym_file)           # bilat.symmetry object

# ---------------------------
# BASIC SANITY CHECKS
# ---------------------------
stopifnot(dim(coords)[3] == nrow(meta))
stopifnot(all(c("specimen_uid","side","series") %in% names(meta)))

cat("Loaded", dim(coords)[3], "wing images\n")

# ============================================================
# PART 1 — CONSENSUS LANDMARK CONFIGURATION
# ============================================================

keep <- meta$series == 1
coords_plot <- coords[,,keep, drop = FALSE]

cat("Consensus plot uses", dim(coords_plot)[3], "wings (series 1 only)\n")

gpa <- gpagen(coords_plot, print.progress = FALSE)
Y   <- gpa$coords
mean_shape <- mshape(Y)

# Flip vertically for display only
mean_shape[,2] <- -mean_shape[,2]
Y[,2,] <- -Y[,2,]

# Plot limits (tight, padded)
pad <- 0.05
xlim <- range(mean_shape[,1]) + c(-pad, pad)
ylim <- range(mean_shape[,2]) + c(-pad, pad)

links <- matrix(c(
  1,2, 2,3, 3,4, 4,5, 5,1,
  3,11, 11,12, 12,13,
  11,10, 10,9, 9,4,
  10,14, 14,15, 15,16,
  16,18,
  18,19,
  19,20,
  20,13,
  20,14,
  9,8, 8,7, 7,6, 6,5,
  8,17, 17,16
), byrow = TRUE, ncol = 2)

png(file.path(out_dir, "05B_consensus_landmarks.png"),
    width = 4000, height = 2400, res = 300)

plot(mean_shape[,1], mean_shape[,2],
     type = "n",
     asp = 1,
     xlim = xlim,
     ylim = ylim,
     xlab = "",
     ylab = "",
     main = "Consensus forewing landmark configuration\n(after Procrustes superimposition)")

# All specimens
for (i in seq_len(dim(Y)[3])) {
  points(Y[,1,i], Y[,2,i], pch = 16, cex = 0.35, col = "grey40")
}

# Mean landmarks
points(mean_shape[,1], mean_shape[,2],
       pch = 21, bg = "white", cex = 1.6)

# Links
for (i in seq_len(nrow(links))) {
  segments(mean_shape[links[i,1],1], mean_shape[links[i,1],2],
           mean_shape[links[i,2],1], mean_shape[links[i,2],2])
}

# Landmark labels
text(mean_shape[,1], mean_shape[,2],
     labels = seq_len(nrow(mean_shape)),
     pos = 3, cex = 1.0)

dev.off()


# ============================================================
# PART 2 — PCA (Symmetric shape and FA shape)
# ============================================================

library(geomorph)

pca_dir <- file.path(out_dir, "PCA")
dir.create(pca_dir, showWarnings = FALSE, recursive = TRUE)

# Helper: keep only specimens with complete landmark data (no NA)
drop_na_specimens <- function(A, ids = NULL) {
  # A is p x k x n
  ok <- apply(A, 3, function(x) all(is.finite(x)))
  A2 <- A[, , ok, drop = FALSE]
  ids2 <- if (!is.null(ids)) ids[ok] else NULL
  list(A = A2, ids = ids2, ok = ok)
}

# Helper: percent variance for axis labels
pc_pct <- function(pca, i) round(100 * (pca$d[i]^2) / sum(pca$d^2), 1)

# --------------------------------------------
# Build specimen-level metadata (unique by uid)
# (important: bilat.symmetry outputs specimen-level arrays)
# --------------------------------------------
meta_uid <- meta[meta$series == 1 & meta$side == "Left",  # any consistent subset
                 c("specimen_uid", "urbanity", "species", "site")]
meta_uid <- meta_uid[!duplicated(meta_uid$specimen_uid), ]

# --------------------------------------------
# 2A) PCA of symmetric component
# sym$symm.shape is p x k x n_spec
# --------------------------------------------
Y_sym <- sym$symm.shape
uid_sym <- dimnames(Y_sym)[[3]]
if (is.null(uid_sym)) uid_sym <- names(sym$unsigned.AI)  # fallback

tmp <- drop_na_specimens(Y_sym, uid_sym)
Y_sym2 <- tmp$A
uid_sym2 <- tmp$ids

stopifnot(dim(Y_sym2)[3] >= 3)  # need at least 3 to PCA meaningfully

pca_sym <- gm.prcomp(Y_sym2)   # correct input type :contentReference[oaicite:5]{index=5}

# Plot using plot.gm.prcomp (recommended; replaces plotTangentSpace) :contentReference[oaicite:6]{index=6}
grp_sym <- factor(meta_uid$urbanity[match(uid_sym2, meta_uid$specimen_uid)],
                  levels = c("Rural", "Urban"))

P <- plot(pca_sym,
          pch = 21,
          bg  = ifelse(grp_sym == "Urban", "grey70", "white"),
          main = "PCA — Symmetric shape")

# Add legend on top of the plot
legend("topright", legend = c("Urban", "Rural"),
       pch = 21, pt.bg = c("grey70", "white"), bty = "n")

# Save a static PNG by replaying the plot call
png(file.path(pca_dir, "05B_PCA_symmetric_shape.png"), width = 1600, height = 1200, res = 200)
plot(pca_sym,
     pch = 21,
     bg  = ifelse(grp_sym == "Urban", "grey70", "white"),
     main = "PCA — Symmetric shape")
legend("topright", legend = c("Urban", "Rural"),
       pch = 21, pt.bg = c("grey70", "white"), bty = "n")
dev.off()

# --------------------------------------------
# 2B) PCA of FA component (asymmetric shape)
# sym$FA.component is p x k x n_spec
# --------------------------------------------
Y_fa <- sym$FA.component
uid_fa <- dimnames(Y_fa)[[3]]
if (is.null(uid_fa)) uid_fa <- names(sym$unsigned.AI)

tmp <- drop_na_specimens(Y_fa, uid_fa)
Y_fa2 <- tmp$A
uid_fa2 <- tmp$ids

stopifnot(dim(Y_fa2)[3] >= 3)

pca_fa <- gm.prcomp(Y_fa2)

grp_fa <- factor(meta_uid$urbanity[match(uid_fa2, meta_uid$specimen_uid)],
                 levels = c("Rural", "Urban"))

png(file.path(pca_dir, "05B_PCA_FA_component.png"), width = 1600, height = 1200, res = 200)
plot(pca_fa,
     pch = 21,
     bg  = ifelse(grp_fa == "Urban", "grey70", "white"),
     main = "PCA — FA (asymmetric) component")
legend("topright", legend = c("Urban", "Rural"),
       pch = 21, pt.bg = c("grey70", "white"), bty = "n")
dev.off()







































