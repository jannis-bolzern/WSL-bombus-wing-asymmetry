# 1 - IMAGE PROCESSING ----------------------------------------------------
#
# This script:
#   - Reads images with standardized filenames
#   - Determines wing side, type, and series
#   - Mirrors images as needed to standardize orientation
#   - Splits forewings and hindwings into separate folders
#
# Input:
#   raw_images/
#
# Output:
#   processed_forewings/
#   processed_hindwings/
#
# All wing images MUST:
#   - be in .jpg format
#   - follow the naming convention below
#
# Example:
#   02-ZHRD-P-LF1.jpg
#
# Structure:
#   <ID>-<CITY><URBANITY><SITE>-<SPECIES>-<SIDE><WING><VIEW>.jpg
#
# CITY:      ZH (Zurich), BS (Basel), BE (Bern)
# URBANITY:  R (Rural), U (Urban)
# SITE:      A–F
# SPECIES:   L (B. lapidarius), P (B. pascuorum)
# SIDE:      L (Left), R (Right)
# WING:      F (Forewing), H (Hindwing)
# VIEW:      1 (Dorsal), 2 (Ventral)
#
# Files not following this convention will fail.

# 1.1 - Libraries ---------------------------------------------------------

if (!requireNamespace("magick", quietly = TRUE)) {
  install.packages("magick")}

library(magick)

# 1.2 - Filename Validation ----------------------------------------------

is_valid_filename <- function(fname) {
  pattern <- "^[0-9]{2}-(ZH|BS|BE)[RU][A-F]-(L|P)-(L|R)(F|H)[12]\\.jpg$"
  grepl(pattern, fname)
}

# 1.3 - Processing Function -----------------------------------------------

process_bumblebee_wings <- function(
    input_dir = "images/raw_images",
    out_fw = "images/processed_forewings",
    out_hw = "images/processed_hindwings"
) {
  
  dir.create(out_fw, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_hw, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(
    input_dir,
    pattern = "\\.jpg$",
    ignore.case = TRUE,
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    stop("No .jpg images found in: ", input_dir)
  }
  
  cat("Found", length(files), "image files\n")
  
  # Validate filenames
  filenames <- basename(files)
  invalid <- filenames[!is_valid_filename(filenames)]
  
  if (length(invalid) > 0) {
    cat("\nERROR: Invalid filename(s) detected:\n")
    for (f in invalid) cat("  -", f, "\n")
    
    stop(
      "\nProcessing stopped.\n",
      "All image files must follow the naming convention:\n",
      "<ID>-<CITY><URBANITY><SITE>-<SPECIES>-<SIDE><WING><VIEW>.jpg\n",
      "Example: 02-ZHRD-P-LF1.jpg"
    )
  }
  
  cat("Filename format check passed \n\n")
  
  results <- data.frame(
    original = character(0),
    processed = character(0),
    stringsAsFactors = FALSE
  )

  # Process each image
  for (img_path in files) {
    
    filename <- basename(img_path)
    
    # Example format: "02-ZHRD-L-LF1.jpg"
    species   <- sub("^.*-(L|P)-[LR][FH][12]\\.jpg$", "\\1", filename)
    side      <- sub("^.*-([LR])[FH][12]\\..*$", "\\1", filename)
    wing_type <- sub("^.*-[LR]([FH])[12]\\..*$", "\\1", filename)
    series    <- sub("^.*-[LR][FH]([12])\\..*$", "\\1", filename)
    series    <- as.integer(series)
    
    species_map <- c(L = "lapidarius",P = "pascuorum")
    
    # Load image
    img <- image_read(img_path)
    mirror <- FALSE
    
    if (side == "L") mirror <- TRUE
      
    # Apply mirroring if required
    if (mirror) {
      cat("Mirroring:", filename, "\n")
      img <- image_flop(img)
    }
    
    # Select output folder
    base_out <- if (wing_type == "F") out_fw else out_hw
    
    # Species-specific subfolder
    out_dir <- file.path(base_out, unname(species_map[species]))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    out_path <- file.path(out_dir, filename)
    
    # Save processed image
    image_write(img, path = out_path)
    
    results <- rbind(results, data.frame(
      original = img_path,
      processed = out_path,
      stringsAsFactors = FALSE
    ))
  }
  
  cat("\nImage processing complete.\n")
  cat("Forewings saved in:", out_fw, "\n")
  cat("Hindwings saved in:", out_hw, "\n")
  
  return(results)
}

# 1.4 - Run Function ------------------------------------------------------

output_summary <- process_bumblebee_wings()
print(head(output_summary))
