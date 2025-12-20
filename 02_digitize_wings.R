# 2 - LANDMARK DIGITIZATION (StereoMorph) ---------------------------------
#
# This script:
#   - Loads processed wing images
#   - Randomizes digitization order
#   - Launches StereoMorph
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

# 2.1 - Libraries ---------------------------------------------------------

if (!requireNamespace("StereoMorph", quietly = TRUE)) {
  install.packages("StereoMorph")}

library(StereoMorph)

# 2.2 - Digitization Setting ----------------------------------------------

# Digitization mode:
# "skip" = skip already digitized images
# "review" = review previously completed images
digitization_mode <- "review"

# Species to digitize:
# e.g. "lapidarius" or "pascuorum"
species <- "lapidarius"

# Wing to digitize:
# "F" = forewing
# "H" = hindwing
digitize_wing <- "F"

# Number of landmarks to digitize:
fw_landmarks <- 21 # forewing landmarks
hw_landmarks <- 6 # hindwing landmarks

# Replicate digitization settings
rep_percent <- 20 # percent of already-digitized images

# 2.3 - Paths -------------------------------------------------------------

if (digitize_wing == "F") {
  image_folder  <- file.path("images", "processed_forewings", species)
  shapes_folder <- file.path("data", "landmark_data_forewings", species)
  lm_file       <- file.path("data", "landmarks_forewings.txt")
  num_landmarks <- fw_landmarks
} else if (digitize_wing == "H") {
  image_folder  <- file.path("images", "processed_hindwings", species)
  shapes_folder <- file.path("data", "landmark_data_hindwings", species)
  lm_file       <- file.path("data", "landmarks_hindwings.txt")
  num_landmarks <- hw_landmarks
} else {
  stop("digitize_wing must be 'F' or 'H'")
}

dir.create(shapes_folder, recursive = TRUE, showWarnings = FALSE)

# 2.4 - Landmark Template -------------------------------------------------

writeLines(paste0("LM", 1:num_landmarks), lm_file)

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
    pattern = "\\.jpg$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(imgs) == 0) {
    warning("No images found in ", image_folder)
    return(list(images = character(0), shapes = character(0)))
  }
  
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

# 4.6 - Function to prepare replicate digitization ------------------------

prepare_replicate_digitization <- function(
    image_folder,
    shapes_folder,
    percent = 20,
    seed = 123,
    mode = "skip"
) {
  
  set.seed(seed)
  
  # Series-1 images only
  imgs1 <- list.files(
    image_folder,
    pattern = "1\\.jpg$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(imgs1) == 0) {
    warning("No series-1 images found for replicate digitization.")
    return(list(images = character(0), shapes = character(0)))
  }
  
  # Series-1 shapes must exist
  shapes1 <- file.path(
    shapes_folder,
    sub("1\\.jpg$", "1.txt", basename(imgs1))
  )
  
  has_dig1 <- file.exists(shapes1)
  imgs1    <- imgs1[has_dig1]
  
  if (length(imgs1) == 0) {
    warning("No existing series-1 digitizations found.")
    return(list(images = character(0), shapes = character(0)))
  }
  
  # Corresponding series-2 shapes
  shapes2 <- file.path(
    shapes_folder,
    sub("1\\.jpg$", "2.txt", basename(imgs1))
  )
  
  exists2 <- file.exists(shapes2)
  
  if (mode == "review") {
    cat("REVIEW mode: showing ONLY existing replicate digitizations\n")
    imgs_sel  <- imgs1[exists2]
    shapes2   <- shapes2[exists2]
  } else {
    # skip mode = create missing replicates
    n_total <- length(imgs1)
    n_rep   <- max(1, ceiling(n_total * percent / 100))
    
    sel <- sample(seq_len(n_total), n_rep)
    
    imgs_sel <- imgs1[sel]
    shapes2  <- shapes2[sel]
    
    # Skip existing series-2
    exists2 <- file.exists(shapes2)
    if (any(exists2)) {
      cat("Skipping", sum(exists2), "already-digitized series-2 files\n")
    }
    
    imgs_sel <- imgs_sel[!exists2]
    shapes2  <- shapes2[!exists2]
  }
  
  if (length(imgs_sel) == 0) {
    warning("No replicate digitizations selected for this mode.")
    return(list(images = character(0), shapes = character(0)))
  }
  
  list(images = imgs_sel, shapes = shapes2)
}


# 2.6 - Run StereoMorph Digitization --------------------------------------

digitize_input <- prepare_digitize_vectors(
  image_folder, shapes_folder, mode = digitization_mode
)

if (length(digitize_input$images) > 0) {
  wing_label <- ifelse(digitize_wing == "F", "FOREWINGS", "HINDWINGS")
  
  cat("\nLaunching StereoMorph for", wing_label, "(",
      length(digitize_input$images), " images )\n")
  
  digitizeImages(
    image.file    = digitize_input$images,
    shapes.file   = digitize_input$shapes,
    landmarks.ref = lm_file
  )
  
  cat(wing_label, "digitization complete.\n")
  invisible(readline())
} else {
  cat(wing_label, "digitization skipped (nothing to do in this mode).\n")
}

# 2.7 - Replicate Digitization --------------------------------------------
  
cat("\nPreparing replicate digitization (", rep_percent, "% of already-digitized wings )\n")
  
rep_input <- prepare_replicate_digitization(
  image_folder, shapes_folder, percent = rep_percent, mode = digitization_mode
)

if (length(rep_input$images) > 0) {
  wing_label <- ifelse(digitize_wing == "F", "FOREWINGS", "HINDWINGS")
    
  cat("Launching StereoMorph for replicate digitization of", wing_label, "(",
      length(rep_input$images), " images )\n")
  
  digitizeImages(
    image.file    = rep_input$images,
    shapes.file   = rep_input$shapes,
    landmarks.ref = lm_file
    )
    
  cat("Replicate digitization complete.\n")
  invisible(readline())
  
  } else {
    cat("No replicate digitizations to perform.\n")
}
