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
#
# Requires:
#   library(StereoMorph)
#   library(geomorph)


# 3.1 - Libraries ---------------------------------------------------------

library(StereoMorph)
library(geomorph)

# 3.2 - Paths -------------------------------------------------------------

fw_shapes <- "landmark_data_forewings"
hw_shapes <- "landmark_data_hindwings"


# 3.3 - Parse Metadata ----------------------------------------------------

# Example filename:
#       02-ZHRD-LF1.txt
#
# Structure:
#   02   = specimen ID
#   ZH   = city
#   R    = rural / U = urban
#   D    = site
#
#   After dash:
#     L  = species (Lapidarius) / P = Pascuorum
#     L/R = left/right
#     F/H = forewing/hindwing
#     1/2 = dorsal/ventral

parse_metadata <- function(fname) {
  
  base <- tools::file_path_sans_ext(basename(fname))
  
  pattern <- "^([0-9]{2})-([A-Z]{2})([RU])([A-F])-([LP])-([LR])([FH])([12])$"
  
  x <- regexec(pattern, base)
  m <- regmatches(base, x)[[1]]
  
  if (length(m) == 0) {
    stop("Filename does not match expected structure: ", fname)
  }
  
  id        <- m[2]
  city      <- m[3]
  urbanity  <- m[4]
  site      <- m[5]
  species   <- m[6]
  side      <- m[7]
  wing      <- m[8]
  view      <- m[9]
  
  urbanity_full <- ifelse(urbanity == "R", "Rural", "Urban")
  
  species_full <- ifelse(
    species == "L", "Bombus lapidarius",
    "Bombus pascuorum"
  )
  
  wing_full <- ifelse(
    wing == "F", "Forewing",
    "Hindwing"
  )
  
  view_full <- ifelse(
    view == "1", "Dorsal",
    "Ventral"
  )
  
  side_full <- ifelse(side == "L", "Left", "Right")
  
  specimen_uid <- paste(id, site, species, sep = "_")
  
  list(
    specimen_uid = specimen_uid,
    specimen_id  = id,
    city         = city,
    urbanity     = urbanity_full,
    site         = paste0(city, "-", site),
    location     = paste0(city, "-", urbanity_full),
    species_code = species,
    species      = species_full,
    side         = side_full,
    wing         = wing_full,
    series       = as.integer(view),
    view         = view_full
  )
}


# 3.4 - Read Shapes -------------------------------------------------------

read_shapes_to_array <- function(folder) {
  
  files <- list.files(folder, pattern = "\\.txt$", full.names = TRUE)
  if (length(files) == 0) stop("No .txt shape files found in: ", folder)
  
  first_shp <- NULL
  for (f in files) {
    shp_i <- tryCatch(readShapes(f), error = function(e) NULL)
    if (!is.null(shp_i)) {
      # probably not needed?
      if (!is.null(shp_i$landmarks.scaled)) {
        first_shp <- shp_i$landmarks.scaled
      } else if (!is.null(shp_i$landmarks.pixel)) {
        first_shp <- shp_i$landmarks.pixel
      }
      if (!is.null(first_shp)) break
    }
  }
  if (is.null(first_shp)) stop("Could not read landmarks from any file in: ", folder)
  
  # ensure first_shp is matrix-like
  first_mat <- as.matrix(first_shp)
  p <- nrow(first_mat)
  k <- ncol(first_mat)
  if (k != 2) stop("Expected 2 columns (x,y) in landmarks, found: ", k)
  
  n <- length(files)
  coords <- array(NA_real_, dim = c(p, k, n))
  
  metadata <- data.frame(
    file = basename(files),
    specimen_uid = character(n),
    specimen_id = character(n),
    city = character(n),
    urbanity = character(n),
    site = character(n),
    location = character(n),
    species_code = character(n),
    species = character(n),
    side = character(n),
    wing = character(n),
    series = integer(n),
    view = character(n),
    stringsAsFactors = FALSE
  )
  
  # iterate files and populate
  for (i in seq_along(files)) {
    f <- files[i]
    shp <- tryCatch(readShapes(f), error = function(e) {
      warning("Failed reading shapes for file: ", basename(f), " — skipping. Error: ", conditionMessage(e))
      return(NULL)
    })
    if (is.null(shp)) next
    
    # pick landmarks: prefer scaled then pixel
    lm_mat <- NULL
    if (!is.null(shp$landmarks.scaled)) {
      lm_mat <- shp$landmarks.scaled
    } else if (!is.null(shp$landmarks.pixel)) {
      lm_mat <- shp$landmarks.pixel
    } else {
      warning("No landmarks found in file: ", basename(f), "; leaving NA for this specimen.")
      next
    }
    
    # coerce to numeric matrix with two columns
    lm_mat <- as.matrix(lm_mat)
    # If the matrix has more than 2 columns, try first two; else error
    if (ncol(lm_mat) < 2) {
      warning("Landmark matrix for ", basename(f), " has <2 columns; skipping.")
      next
    } else if (ncol(lm_mat) > 2) {
      lm_mat <- lm_mat[, 1:2]
    }
    
    # If landmark count differs from 'p', pad or truncate gracefully
    p_this <- nrow(lm_mat)
    if (p_this < p) {
      # create p x 2 with NA, then fill top p_this rows
      tmp <- matrix(NA_real_, nrow = p, ncol = 2)
      tmp[seq_len(p_this), ] <- apply(lm_mat, 2, as.numeric)
      lm_mat_num <- tmp
    } else if (p_this > p) {
      # truncate to first p rows (warn)
      warning("File ", basename(f), " has ", p_this, " landmarks (expected ", p, "). Truncating to first ", p, ".")
      lm_mat_num <- apply(lm_mat[seq_len(p), , drop = FALSE], 2, as.numeric)
    } else {
      lm_mat_num <- apply(lm_mat, 2, as.numeric)
    }
    
    coords[, , i] <- lm_mat_num
    
    # parse metadata using your parse_metadata() function (must exist in environment)
    md <- tryCatch(parse_metadata(basename(f)), error = function(e) {
      warning("parse_metadata() failed for ", basename(f), ": ", conditionMessage(e))
      return(NULL)
    })
    if (!is.null(md)) {
      metadata$specimen_uid[i] <- md$specimen_uid 
      metadata$specimen_id[i]  <- md$specimen_id
      metadata$city[i]         <- md$city
      metadata$urbanity[i]     <- md$urbanity
      metadata$site[i]         <- md$site
      metadata$location[i]     <- md$location
      metadata$species_code[i] <- md$species_code
      metadata$species[i]      <- md$species
      metadata$side[i]         <- md$side
      metadata$wing[i]         <- md$wing
      metadata$series[i]       <- md$series
      metadata$view[i]         <- md$view
    }
  }
  
  return(list(coords = coords, metadata = metadata))
}

# 3.5 - Process Forewings -------------------------------------------------

cat("Processing FOREWING shape files...\n")
fw <- read_shapes_to_array(fw_shapes)

saveRDS(fw$coords, file = "fw_coords.RDS")
write.csv(fw$metadata, file = "fw_metadata.csv", row.names = FALSE)

cat("Forewing data saved:\n  fw_coords.RDS\n  fw_metadata.csv\n\n")

# 3.6 - Process Hindwings -------------------------------------------------

cat("Processing HINDWING shape files...\n")
hw <- read_shapes_to_array(hw_shapes)

saveRDS(hw$coords, file = "hw_coords.RDS")
write.csv(hw$metadata, file = "hw_metadata.csv", row.names = FALSE)

cat("Hindwing data saved:\n  hw_coords.RDS\n  hw_metadata.csv\n\n")

# 3.6 - Summary -----------------------------------------------------------

cat("Done!\n")
cat("Forewings:", dim(fw$coords)[3], "specimens\n")
cat("Hindwings:", dim(hw$coords)[3], "specimens\n")
