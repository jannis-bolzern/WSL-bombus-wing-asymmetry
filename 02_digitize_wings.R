# ============================================================
# 02 — DIGITIZATION SCRIPT (StereoMorph)
# ============================================================
# This script:
#   - Loads processed wing images
#   - Randomizes digitization order
#   - Skips already-digitized wings
#   - Launches StereoMorph for forewings and hindwings separately
#
# Input:
#   processed_forewings/
#   processed_hindwings/
#
# Output:
#   landmark_data_forewings/
#   landmark_data_hindwings/
#
# Requires:
#   install.packages("StereoMorph")
# ============================================================

library(StereoMorph)

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------
fw_folder   <- "processed_forewings"
hw_folder   <- "processed_hindwings"

fw_shapes   <- "landmark_data_forewings"
hw_shapes   <- "landmark_data_hindwings"

if (!dir.exists(fw_shapes)) dir.create(fw_shapes)
if (!dir.exists(hw_shapes)) dir.create(hw_shapes)

# ------------------------------------------------------------
# LANDMARK TEMPLATES
# ------------------------------------------------------------
fw_num_landmarks <- 20
fw_lm_file <- "landmarks_forewings.txt"
writeLines(paste0("LM", 1:fw_num_landmarks), fw_lm_file)

hw_num_landmarks <- 6
hw_lm_file <- "landmarks_hindwings.txt"
writeLines(paste0("LM", 1:hw_num_landmarks), hw_lm_file)

# ------------------------------------------------------------
# Helper: prepare digitization lists
# ------------------------------------------------------------
prepare_digitize_vectors <- function(
    image_folder, 
    shapes_folder, 
    ext_out = ".txt",
    randomize = TRUE,
    overwrite = FALSE
) {
  imgs <- list.files(
    image_folder,
    pattern = "\\.(jpg|jpeg|png|tif|tiff)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(imgs) == 0)
    stop("No images found in ", image_folder)
  
  if (randomize)
    imgs <- sample(imgs)  # randomize order
  
  shapes_paths <- file.path(
    shapes_folder,
    paste0(tools::file_path_sans_ext(basename(imgs)), ext_out)
  )
  
  # Skip images already digitized
  if (!overwrite) {
    exists <- file.exists(shapes_paths)
    if (any(exists)) {
      cat("Skipping", sum(exists), "already-digitized images:\n")
      print(basename(imgs[exists]))
      
      imgs <- imgs[!exists]
      shapes_paths <- shapes_paths[!exists]
      
      if (length(imgs) == 0) {
        warning("All images already digitized.")
        return(list(images = character(0), shapes = character(0)))
      }
    }
  }
  
  return(list(images = imgs, shapes = shapes_paths))
}

# ------------------------------------------------------------
# FOREWING DIGITIZATION
# ------------------------------------------------------------
fw_input <- prepare_digitize_vectors(fw_folder, fw_shapes)

if (length(fw_input$images) > 0) {
  cat("\nLaunching StereoMorph for FOREWINGS (",
      length(fw_input$images), " images )\n")
  
  digitizeImage(
    image.file    = fw_input$images,
    shapes.file   = fw_input$shapes,
    landmarks.ref = fw_lm_file
  )
  
  cat("Forewing digitization complete. Press ENTER to continue.\n")
  invisible(readline())
} else {
  cat("Forewing digitization skipped (nothing to do).\n")
}

# ------------------------------------------------------------
# HINDWING DIGITIZATION
# ------------------------------------------------------------
hw_input <- prepare_digitize_vectors(hw_folder, hw_shapes)

if (length(hw_input$images) > 0) {
  cat("\nLaunching StereoMorph for HINDWINGS (",
      length(hw_input$images), " images )\n")
  
  digitizeImage(
    image.file    = hw_input$images,
    shapes.file   = hw_input$shapes,
    landmarks.ref = hw_lm_file
  )
  
  cat("Hindwing digitization complete.\n")
} else {
  cat("Hindwing digitization skipped (nothing to do).\n")
}
