# 4 - BUILD GEOMORPH COORDINATE ARRAYS ------------------------------------
#
# This script:
#   - Reads StereoMorph .txt landmarks from species subfolders
#   - Builds geomorph arrays (p × k × n) per species
#   - Extracts full metadata from filenames
#   - Saves per-species outputs:
#         data/fw_coords_<SPECIES>.RDS
#         data/fw_metadata_<SPECIES>.csv
#         data/hw_coords_<SPECIES>.RDS
#         data/hw_metadata_<SPECIES>.csv

# 4.1 - Libraries ---------------------------------------------------------

if (!requireNamespace("StereoMorph", quietly = TRUE)) install.packages("StereoMorph")
library(StereoMorph)

if (!requireNamespace("geomorph", quietly = TRUE)) install.packages("geomorph")
library(geomorph)

# 4.2 - Base paths --------------------------------------------------------

fw_base <- "data/landmark_data_forewings"
hw_base <- "data/landmark_data_hindwings"

out_base <- "data"

# 4.3 - Parse Metadata ----------------------------------------------------

parse_metadata <- function(fname) {
  
  base <- tools::file_path_sans_ext(basename(fname))
  pattern <- "^([0-9]{2})-([A-Z]{2})([RU])([A-F])-([LP])-([LR])([FH])([12])$"
  
  x <- regexec(pattern, base)
  m <- regmatches(base, x)[[1]]
  
  if (length(m) == 0) {
    stop("Filename does not match expected structure: ", fname)
  }
  
  id            <- m[2]
  city          <- m[3]
  urbanity_code <- m[4]
  site_letter   <- m[5]
  species_code  <- m[6]
  side_code     <- m[7]
  wing_code     <- m[8]
  series_code   <- m[9]
  
  urbanity_full <- ifelse(urbanity_code == "R", "Rural", "Urban")
  
  species_full <- ifelse(
    species_code == "L", "Bombus lapidarius",
    "Bombus pascuorum"
  )
  
  wing_full <- ifelse(
    wing_code == "F", "Forewing",
    "Hindwing"
  )
  
  side_full <- ifelse(side_code == "L", "Left", "Right")
  
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
    series         = as.integer(series_code)   # digitization replicate
  )
}

# 4.4 - Read Shapes -------------------------------------------------------

read_shapes_to_array <- function(folder) {
  
  if (!dir.exists(folder)) {
    message("  → Skipping (folder does not exist): ", folder)
    return(NULL)
  }
  
  files <- list.files(
    folder,
    pattern = "\\.txt$",
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    message("  → Skipping (no landmark files found): ", folder)
    return(NULL)
  }
  
  base_names <- basename(files)
  if (any(duplicated(base_names))) {
    dups <- base_names[duplicated(base_names)]
    stop(
      "Duplicate .txt filenames found in folder: ", folder,
      "\nExamples: ", paste(head(dups, 10), collapse = ", ")
    )
  }
  
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
    
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(files)) {
    
    shp <- readShapes(files[i])
    lm  <- as.matrix(shp$landmarks.pixel)
    
    coords[, , i] <- lm
    
    md <- parse_metadata(basename(files[i]))
    
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
  }
  
  key <- paste(metadata$specimen_uid, metadata$side, metadata$series, sep = "__")
  if (any(duplicated(key))) {
    stop("Duplicate specimen_uid × side × series in: ", folder)
  }
  
  list(coords = coords, metadata = metadata)
}

# 4.5 - Process per species ----------------------------------------------

fw_species_dirs <- list.dirs(fw_base, recursive = FALSE, full.names = TRUE)
hw_species_dirs <- list.dirs(hw_base, recursive = FALSE, full.names = TRUE)

species_dirs <- sort(unique(c(fw_species_dirs, hw_species_dirs)))
species <- sort(unique(basename(species_dirs)))

for (sp in species) {
  
  cat("\nProcessing species:", sp, "\n")
  
  fw <- read_shapes_to_array(file.path(fw_base, sp))
  if (!is.null(fw)) {
    saveRDS(
      fw$coords,
      file = file.path(out_base, paste0("fw_coords_", sp, ".RDS"))
    )
    write.csv(
      fw$metadata,
      file = file.path(out_base, paste0("fw_metadata_", sp, ".csv")),
      row.names = FALSE
    )
    cat("  Forewings:", dim(fw$coords)[3], "images\n")
  }
  
  hw <- read_shapes_to_array(file.path(hw_base, sp))
  if (!is.null(hw)) {
    saveRDS(
      hw$coords,
      file = file.path(out_base, paste0("hw_coords_", sp, ".RDS"))
    )
    write.csv(
      hw$metadata,
      file = file.path(out_base, paste0("hw_metadata_", sp, ".csv")),
      row.names = FALSE
    )
    cat("  Hindwings:", dim(hw$coords)[3], "images\n")
  }
  cat("Geomorph array construction complete.\n\n")
  
  cat("Outputs written to working directory:\n")
  cat("  - Forewing coordinates:  fw_coords_<SPECIES>.RDS\n")
  cat("  - Forewing metadata:     fw_metadata_<SPECIES>.csv\n")
  cat("  - Hindwing coordinates:  hw_coords_<SPECIES>.RDS\n")
  cat("  - Hindwing metadata:     hw_metadata_<SPECIES>.csv\n\n")
  
  cat("Species processed:\n")
  for (sp in species) {
    if (file.exists(paste0("fw_coords_", sp, ".RDS"))) {
      cat("  -", sp, "\n")
    }
  }
}
