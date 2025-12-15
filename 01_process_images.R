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

# 1.1 - Libraries ---------------------------------------------------------

if (!requireNamespace("magick", quietly = TRUE)) {
  install.packages("magick")}

library(magick)

# 1.2 - Processing Function -----------------------------------------------

process_bumblebee_wings <- function(
    input_dir = "raw_images",
    out_fw = "processed_forewings",
    out_hw = "processed_hindwings"
) {
  
  # Create directories if not already existent
  dir.create(out_fw, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_hw, recursive = TRUE, showWarnings = FALSE)
  
  # Load image file list
  files <- list.files(
    input_dir,
    pattern = "\\.(jpg|jpeg|png|tif|tiff)$",
    ignore.case = TRUE,
    full.names = TRUE
  )
  
  cat("Found", length(files), "image files\n\n")
  
  results <- data.frame(
    original = character(0),
    processed = character(0),
    stringsAsFactors = FALSE
  )
  
  # Process each image
  for (img_path in files) {
    
    filename <- basename(img_path)
    
    # Example format: "02-ZHRD-L-LF1.jpg"
    side      <- sub("^.*-([LR])[FH][12]\\..*$", "\\1", filename)
    wing_type <- sub("^.*-[LR]([FH])[12]\\..*$", "\\1", filename)
    series    <- sub("^.*-[LR][FH]([12])\\..*$", "\\1", filename)
    series    <- as.integer(series)
    
    # Load image
    img <- image_read(img_path)
    mirror <- FALSE
    
    # Series 1 - mirror LEFT wings
    if (series == 1 && side == "L") mirror <- TRUE
    
    # Series 2 - mirror RIGHT wings (inverted view)
    if (series == 2 && side == "R") mirror <- TRUE
    
    # Apply mirroring if required
    if (mirror) {
      cat("Mirroring:", filename, "\n")
      img <- image_flop(img)
    }
    
    # Select output folder
    out_dir <- if (wing_type == "F") out_fw else out_hw
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

# 1.3 - Run Function ------------------------------------------------------

output_summary <- process_bumblebee_wings()
print(head(output_summary))
