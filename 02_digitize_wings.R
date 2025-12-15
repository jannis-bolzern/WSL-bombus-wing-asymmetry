# 2 - LANDMARK DIGITIZATION (StereoMorph) ---------------------------------
#
# This script:
#   - Loads processed wing images
#   - Randomizes digitization order
#   - Launches StereoMorph for forewings and hindwings separately
#   - Has 2 modes:
#         "skip"      = skip images already digitized (default)
#         "review"    = ONLY show already-digitized images (QC)
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


# 2.1 - Libraries ---------------------------------------------------------

library(StereoMorph)

# 2.2 - Digitization Setting ----------------------------------------------

# "skip"      → skip already digitized images
# "review"    → ONLY digitize previously completed images

digitization_mode <- "review"

cat("\nDigitization mode set to:", toupper(digitization_mode), "\n\n")

# 2.3 - Paths -------------------------------------------------------------

fw_folder   <- "processed_forewings"
hw_folder   <- "processed_hindwings"

fw_shapes   <- "landmark_data_forewings"
hw_shapes   <- "landmark_data_hindwings"

if (!dir.exists(fw_shapes)) dir.create(fw_shapes)
if (!dir.exists(hw_shapes)) dir.create(hw_shapes)

# 2.4 - Landmark Templates ------------------------------------------------

fw_num_landmarks <- 21
fw_lm_file <- "landmarks_forewings.txt"
writeLines(paste0("LM", 1:fw_num_landmarks), fw_lm_file)

hw_num_landmarks <- 6
hw_lm_file <- "landmarks_hindwings.txt"
writeLines(paste0("LM", 1:hw_num_landmarks), hw_lm_file)

# 2.5 - Prepare Digitization Lists ----------------------------------------

prepare_digitize_vectors <- function(
    image_folder, 
    shapes_folder, 
    ext_out = ".txt",
    randomize = TRUE,
    mode = "skip"
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
    imgs <- sample(imgs)
  
  shapes_paths <- file.path(
    shapes_folder,
    paste0(tools::file_path_sans_ext(basename(imgs)), ext_out)
  )
  
  exists <- file.exists(shapes_paths)

  if (mode == "skip") {
    if (any(exists)) {
      cat("Skipping", sum(exists), "already-digitized images:\n")
      print(basename(imgs[exists]))
    }
    imgs <- imgs[!exists]
    shapes_paths <- shapes_paths[!exists]
  }
  
  if (mode == "review") {
    cat("REVIEW mode: showing ONLY already-digitized images\n")
    imgs <- imgs[exists]
    shapes_paths <- shapes_paths[exists]
  }
  
  if (length(imgs) == 0) {
    warning("No images selected for this mode.")
    return(list(images = character(0), shapes = character(0)))
  }
  
  return(list(images = imgs, shapes = shapes_paths))
}

# 2.6 - Run Forewing Digitization -----------------------------------------

fw_input <- prepare_digitize_vectors(
  fw_folder, fw_shapes, mode = digitization_mode
)

if (length(fw_input$images) > 0) {
  cat("\nLaunching StereoMorph for FOREWINGS (",
      length(fw_input$images), " images )\n")
  
  digitizeImages(
    image.file    = fw_input$images,
    shapes.file   = fw_input$shapes,
    landmarks.ref = fw_lm_file
  )
  
  cat("Forewing digitization complete.\n")
  invisible(readline())
} else {
  cat("Forewing digitization skipped (nothing to do in this mode).\n")
}

# 2.7 - Run Hindwing Digitization -----------------------------------------

hw_input <- prepare_digitize_vectors(
  hw_folder, hw_shapes, mode = digitization_mode
)

if (length(hw_input$images) > 0) {
  cat("\nLaunching StereoMorph for HINDWINGS (",
      length(hw_input$images), " images )\n")
  
  digitizeImages(
    image.file    = hw_input$images,
    shapes.file   = hw_input$shapes,
    landmarks.ref = hw_lm_file
  )
  
  cat("Hindwing digitization complete.\n")
} else {
  cat("Hindwing digitization skipped (nothing to do in this mode).\n")
}

