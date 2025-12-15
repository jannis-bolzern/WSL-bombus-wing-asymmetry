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

if (!requireNamespace("StereoMorph", quietly = TRUE)) {
  install.packages("StereoMorph")}

library(StereoMorph)

if (!requireNamespace("geomorph", quietly = TRUE)) {
  install.packages("geomorph")}

library(geomorph)

# 3.2 - Paths -------------------------------------------------------------

fw_shapes <- "landmark_data_forewings"
hw_shapes <- "landmark_data_hindwings"

# 3.3 - Parse Metadata ----------------------------------------------------

# Example filename:
#   02-ZHRD-P-LF1.jpg
#
# Structure:
#   <ID>-<CITY><URBANITY><SITE>-<SPECIES>-<SIDE><WING><VIEW>.jpg
#
# CITY:      ZH (Zurich), BS (Basel), GE (Geneva)
# URBANITY:  R (Rural), U (Urban)
# SITE:      A–F
# SPECIES:   L (B. lapidarius), P (B. pascuorum)
# SIDE:      L (Left), R (Right)
# WING:      F (Forewing), H (Hindwing)
# VIEW:      1 (Dorsal), 2 (Ventral)
#
# - The specimen number alone is NOT a unique identifier.
# - Uniqueness is ensured later via a constructed specimen_uid.
# - Filenames not matching this pattern will be rejected.

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
  if (length(files) == 0)
    stop("No .txt shape files found in: ", folder)
  
  # Read first file to get dimensions
  shp0 <- readShapes(files[1])
  
  if (is.null(shp0$landmarks.pixel))
    stop("No pixel landmarks found in first file: ", basename(files[1]))
  
  first_mat <- as.matrix(shp0$landmarks.pixel)
  
  if (ncol(first_mat) != 2)
    stop("Expected 2 columns (x, y) in landmarks.")
  
  p <- nrow(first_mat)
  k <- 2
  n <- length(files)
  
  coords <- array(NA_real_, dim = c(p, k, n))
  
  metadata <- data.frame(
    file = basename(files),
    specimen_uid = character(n),
    specimen_id  = character(n),
    city         = character(n),
    urbanity     = character(n),
    site         = character(n),
    location     = character(n),
    species_code = character(n),
    species      = character(n),
    side         = character(n),
    wing         = character(n),
    series       = integer(n),
    view         = character(n),
    stringsAsFactors = FALSE
  )
  
  # Loop through files
  for (i in seq_along(files)) {
    f <- files[i]
    
    shp <- tryCatch(
      readShapes(f),
      error = function(e) stop("Failed reading: ", basename(f))
    )
    
    lm <- shp$landmarks.pixel
    if (is.null(lm))
      stop("No pixel landmarks in file: ", basename(f))
    
    lm <- as.matrix(lm)
    
    if (nrow(lm) != p || ncol(lm) != 2)
      stop("Landmark dimensions differ in file: ", basename(f))
    
    coords[, , i] <- lm
    
    md <- parse_metadata(basename(f))
    
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

# 3.6 - Summary -----------------------------------------------------------

cat("Forewings:", dim(fw$coords)[3], "specimens\n")
cat("Hindwings:", dim(hw$coords)[3], "specimens\n")
