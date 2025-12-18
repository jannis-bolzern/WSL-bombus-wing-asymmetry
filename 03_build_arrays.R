# 3 - BUILD GEOMORPH COORDINATE ARRAYS ------------------------------------
#
# This script:
#   - Reads StereoMorph .txt landmarks from forewings & hindwings
#   - Builds geomorph arrays (p × k × n)
#   - Extracts full metadata from filenames
#   - Saves:
#         fw_coords.RDS
#         hw_coords.RDS
#         fw_metadata.csv
#         hw_metadata.csv

# 3.1 - Libraries ---------------------------------------------------------

if (!requireNamespace("StereoMorph", quietly = TRUE)) install.packages("StereoMorph")
library(StereoMorph)

if (!requireNamespace("geomorph", quietly = TRUE)) install.packages("geomorph")
library(geomorph)

# 3.2 - Paths -------------------------------------------------------------

fw_shapes <- "landmark_data_forewings"
hw_shapes <- "landmark_data_hindwings"

# 3.3 - Parse Metadata ----------------------------------------------------
#
# Example filename:
#   02-ZHRD-P-LF1.txt   (same base as jpg; extension doesn't matter)
#
# Structure:
#   <ID>-<CITY><URBANITY><SITE>-<SPECIES>-<SIDE><WING><VIEW>
#
# ID:        01–99
# CITY:      ZH, BS, BE
# URBANITY:  R, U
# SITE:      A–F
# SPECIES:   L, P
# SIDE:      L, R
# WING:      F, H
# VIEW:      1, 2
#
# Note:
#   specimen_uid MUST be unique across all cities/sites/species.
#   We enforce this by constructing:
#     specimen_uid = CITY_URBANITY_SITELETTER_ID_SPECIES
#
# This makes Script 4/5 safe, because pairing uses specimen_uid × side × series.

parse_metadata <- function(fname) {
  
  base <- tools::file_path_sans_ext(basename(fname))
  pattern <- "^([0-9]{2})-([A-Z]{2})([RU])([A-F])-([LP])-([LR])([FH])([12])$"
  
  x <- regexec(pattern, base)
  m <- regmatches(base, x)[[1]]
  
  if (length(m) == 0) {
    stop("Filename does not match expected structure: ", fname)
  }
  
  id            <- m[2]  # 01–99
  city          <- m[3]  # ZH / BS / BE
  urbanity_code <- m[4]  # R / U
  site_letter   <- m[5]  # A-F
  species_code  <- m[6]  # L / P
  side_code     <- m[7]  # L / R
  wing_code     <- m[8]  # F / H
  view_code     <- m[9]  # 1 / 2
  
  urbanity_full <- ifelse(urbanity_code == "R", "Rural", "Urban")
  
  species_full <- ifelse(
    species_code == "L", "Bombus lapidarius",
    "Bombus pascuorum"
  )
  
  wing_full <- ifelse(
    wing_code == "F", "Forewing",
    "Hindwing"
  )
  
  view_full <- ifelse(
    view_code == "1", "Dorsal",
    "Ventral"
  )
  
  side_full <- ifelse(side_code == "L", "Left", "Right")
  
  # Robust unique ID (safe across cities and sites)
  specimen_uid <- paste(city, urbanity_code, site_letter, id, species_code, sep = "_")
  
  list(
    specimen_uid   = specimen_uid,
    specimen_id    = id,
    city           = city,
    urbanity_code  = urbanity_code,
    urbanity       = urbanity_full,
    site_letter    = site_letter,
    site           = paste0(city, "-", site_letter),
    location       = paste0(city, "-", urbanity_full),
    species_code   = species_code,
    species        = species_full,
    side           = side_full,
    wing           = wing_full,
    series         = as.integer(view_code), # replicate index 1/2
    view           = view_full
  )
}

# 3.4 - Read Shapes -------------------------------------------------------

read_shapes_to_array <- function(folder) {
  
  files <- list.files(folder, pattern = "\\.txt$", full.names = TRUE)
  if (length(files) == 0) stop("No .txt shape files found in: ", folder)
  
  # Ensure filenames are unique (basic hygiene)
  base_names <- basename(files)
  if (any(duplicated(base_names))) {
    dups <- base_names[duplicated(base_names)]
    stop("Duplicate .txt filenames found in folder: ", folder, "\nExamples: ", paste(head(dups, 10), collapse = ", "))
  }
  
  # Read first file to get landmark dimensions
  shp0 <- readShapes(files[1])
  if (is.null(shp0$landmarks.pixel)) stop("No pixel landmarks found in first file: ", basename(files[1]))
  
  first_mat <- as.matrix(shp0$landmarks.pixel)
  if (ncol(first_mat) != 2) stop("Expected 2 columns (x, y) in landmarks.")
  
  p <- nrow(first_mat)
  k <- 2
  n <- length(files)
  
  coords <- array(NA_real_, dim = c(p, k, n))
  
  metadata <- data.frame(
    file = basename(files),
    
    specimen_uid  = character(n),
    specimen_id   = character(n),
    city          = character(n),
    urbanity_code = character(n),
    urbanity      = character(n),
    site_letter   = character(n),
    site          = character(n),
    location      = character(n),
    
    species_code  = character(n),
    species       = character(n),
    side          = character(n),
    wing          = character(n),
    series        = integer(n),
    view          = character(n),
    
    stringsAsFactors = FALSE
  )
  
  # Loop through files
  for (i in seq_along(files)) {
    f <- files[i]
    
    shp <- tryCatch(
      readShapes(f),
      error = function(e) stop("Failed reading: ", basename(f), "\n", e$message)
    )
    
    lm <- shp$landmarks.pixel
    if (is.null(lm)) stop("No pixel landmarks in file: ", basename(f))
    
    lm <- as.matrix(lm)
    if (nrow(lm) != p || ncol(lm) != 2) {
      stop("Landmark dimensions differ in file: ", basename(f),
           " (expected ", p, "x2, got ", nrow(lm), "x", ncol(lm), ")")
    }
    
    coords[, , i] <- lm
    
    md <- parse_metadata(basename(f))
    
    metadata$specimen_uid[i]  <- md$specimen_uid
    metadata$specimen_id[i]   <- md$specimen_id
    metadata$city[i]          <- md$city
    metadata$urbanity_code[i] <- md$urbanity_code
    metadata$urbanity[i]      <- md$urbanity
    metadata$site_letter[i]   <- md$site_letter
    metadata$site[i]          <- md$site
    metadata$location[i]      <- md$location
    metadata$species_code[i]  <- md$species_code
    metadata$species[i]       <- md$species
    metadata$side[i]          <- md$side
    metadata$wing[i]          <- md$wing
    metadata$series[i]        <- md$series
    metadata$view[i]          <- md$view
  }
  
  # --- sanity checks -----------------------------------------------------
  
  # 1) specimen_uid must exist
  if (anyNA(metadata$specimen_uid) || any(metadata$specimen_uid == "")) {
    stop("specimen_uid parsing failed for at least one file in: ", folder)
  }
  
  # 2) No duplicates of specimen_uid × side × series (critical for pairing)
  key <- paste(metadata$specimen_uid, metadata$side, metadata$series, sep = "__")
  dup <- key[duplicated(key)]
  if (length(dup) > 0) {
    bad_files <- metadata$file[key %in% dup]
    stop(
      "Duplicate specimen_uid × side × series detected in: ", folder, "\n",
      "This usually means duplicate files or filename collisions.\n",
      "Examples of duplicated keys: ", paste(head(unique(dup), 10), collapse = ", "), "\n",
      "Example files: ", paste(head(bad_files, 10), collapse = ", ")
    )
  }
  
  # 3) Informative warning: missing expected combinations (not an error here)
  #    Script 4 will filter to complete specimens anyway.
  tab <- with(metadata, table(specimen_uid, side, series))
  incomplete <- sum(apply(tab, 1, function(x) any(x == 0)))
  if (incomplete > 0) {
    message("Note: ", incomplete, " specimen_uid entries are missing at least one side/series in ", folder,
            " (this is OK; Script 4 filters to complete L/R×1/2).")
  }
  
  list(coords = coords, metadata = metadata)
}

# 3.5 - Process Forewings -------------------------------------------------

fw <- read_shapes_to_array(fw_shapes)

saveRDS(fw$coords, file = "fw_coords.RDS")
write.csv(fw$metadata, file = "fw_metadata.csv", row.names = FALSE)

# 3.6 - Process Hindwings -------------------------------------------------

hw <- read_shapes_to_array(hw_shapes)

saveRDS(hw$coords, file = "hw_coords.RDS")
write.csv(hw$metadata, file = "hw_metadata.csv", row.names = FALSE)

# 3.7 - Summary -----------------------------------------------------------

cat("Forewings:", dim(fw$coords)[3], "images\n")
cat("Hindwings:", dim(hw$coords)[3], "images\n")
cat("Done.\n")
