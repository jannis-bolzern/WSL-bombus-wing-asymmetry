# 2 - LANDMARK DIGITIZATION (StereoMorph) ---------------------------------
#
# This script:
#   - Loads processed wing images
#   - Creates temporary copies in the image folder for digitization replicates
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
digitization_mode <- "skip"

# Species to digitize:
# e.g. "lapidarius", "pascuorum" or both as c("lapidarius", "pascuorum")
species <- c("lapidarius", "pascuorum")

# Wing to digitize:
# "F" = forewing
# "H" = hindwing
digitize_wing <- "F"

# Number of landmarks to digitize:
fw_landmarks <- 15 # forewing landmarks
hw_landmarks <- 6  # hindwing landmarks

# 2.3 - Paths -------------------------------------------------------------

if (digitize_wing == "F") {
  image_folders  <- file.path("images", "processed_forewings", species)
  shapes_folders <- file.path("data", "landmark_data_forewings", species)
  lm_file       <- file.path("data", "landmarks_forewings.txt")
  num_landmarks <- fw_landmarks
} else if (digitize_wing == "H") {
  image_folders  <- file.path("images", "processed_hindwings", species)
  shapes_folders <- file.path("data", "landmark_data_hindwings", species)
  lm_file       <- file.path("data", "landmarks_hindwings.txt")
  num_landmarks <- hw_landmarks
} else {
  stop("digitize_wing must be 'F' or 'H'")
}

for (sf in shapes_folders) dir.create(sf, recursive = TRUE, showWarnings = FALSE)

# 2.4 - Landmark Template -------------------------------------------------

writeLines(paste0("LM", 1:num_landmarks), lm_file)

# 2.5 – Helper Functions --------------------------------------------------

# Function to clean up temporary files
cleanup_temp_files <- function(temp_files) {
  if (length(temp_files) > 0) {
    unlink(temp_files, recursive = TRUE)
    cat("Cleaned up", length(temp_files), "temporary files\n")
  }
}

# Prepare digitization list for single session
prepare_single_session_jobs <- function(image_folder, shapes_folder, mode = "skip") {
  
  # Initialize vectors
  imgs_to_digitize <- character(0)
  shapes_to_save <- character(0)
  temp_files <- character(0)  # To track temporary files for cleanup
  
  # Get all series 2 images first to know which have imaging replicates
  series2_imgs <- list.files(image_folder, pattern = "2\\.jpg$", full.names = TRUE)
  
  # Create a lookup of which specimens have imaging replicates
  has_imaging_replicate <- character()
  for (img2 in series2_imgs) {
    img1 <- sub("2\\.jpg$", "1.jpg", img2)
    if (file.exists(img1)) {
      has_imaging_replicate <- c(has_imaging_replicate, img1)
    }
  }
  
  # Create a temporary subdirectory in the image folder for replicates
  temp_dir <- file.path(image_folder, "temp_replicates")
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get all series 1 images (original photos)
  series1_imgs <- list.files(image_folder, pattern = "1\\.jpg$", full.names = TRUE)
  
  # Process each series 1 image
  for (img1 in series1_imgs) {
    basename_no_ext <- tools::file_path_sans_ext(basename(img1))
    
    # Regular digitization (image_1.jpg -> image_1.txt)
    shape_regular <- file.path(shapes_folder, paste0(basename_no_ext, ".txt"))
    
    # Check if we should include regular digitization
    if (mode == "skip") {
      include_regular <- !file.exists(shape_regular)
    } else { # mode == "review"
      include_regular <- file.exists(shape_regular)
    }
    
    # Add regular digitization if needed
    if (include_regular) {
      imgs_to_digitize <- c(imgs_to_digitize, img1)
      shapes_to_save <- c(shapes_to_save, shape_regular)
    }
    
    # Only create digitization replicate IF there's a corresponding series 2 image
    if (img1 %in% has_imaging_replicate) {
      shape_rep <- file.path(shapes_folder, paste0(basename_no_ext, "_rep.txt"))
      
      # Check if we should include digitization replicate
      if (mode == "skip") {
        include_rep <- !file.exists(shape_rep)
      } else { # mode == "review"
        include_rep <- file.exists(shape_rep)
      }
      
      # Add digitization replicate if needed
      if (include_rep) {
        # Create temporary copy in the temp directory (same folder as images)
        temp_img <- file.path(temp_dir, paste0(basename_no_ext, "_rep_temp.jpg"))
        
        # Copy the original image to temporary location
        file.copy(img1, temp_img, overwrite = TRUE)
        
        # Add to job list
        imgs_to_digitize <- c(imgs_to_digitize, temp_img)
        shapes_to_save <- c(shapes_to_save, shape_rep)
        
        # Track for cleanup
        temp_files <- c(temp_files, temp_img)
      }
    }
  }
  
  # Process each series 2 image (imaging replicates)
  for (img2 in series2_imgs) {
    # Check if corresponding series 1 exists (only process if part of a pair)
    img1 <- sub("2\\.jpg$", "1.jpg", img2)
    if (!file.exists(img1)) next
    
    basename_no_ext <- tools::file_path_sans_ext(basename(img2))
    shape_regular <- file.path(shapes_folder, paste0(basename_no_ext, ".txt"))
    
    if (mode == "skip") {
      include <- !file.exists(shape_regular)
    } else {
      include <- file.exists(shape_regular)
    }
    
    if (include) {
      imgs_to_digitize <- c(imgs_to_digitize, img2)
      shapes_to_save <- c(shapes_to_save, shape_regular)
    }
  }
  
  # Return list with jobs and temp files for cleanup
  return(list(
    jobs = data.frame(
      image = imgs_to_digitize,
      shape = shapes_to_save,
      stringsAsFactors = FALSE
    ),
    temp_files = temp_files,
    temp_dir = temp_dir
  ))
}

# 2.6 - Run Digitization ---------------------------------------

cat("Digitization mode:", toupper(digitization_mode), "\n")
wing_label <- ifelse(digitize_wing == "F", "Forewings", "Hindwings")
cat("Species:", paste(species, collapse = " + "), "\n")
cat("Wing type:", wing_label, "\n")

# collect jobs from all species folders
all_jobs <- data.frame(image = character(0), shape = character(0), stringsAsFactors = FALSE)
all_temp_files <- character(0)
all_temp_dirs  <- character(0)

for (i in seq_along(species)) {
  prep_result <- prepare_single_session_jobs(
    image_folder  = image_folders[i],
    shapes_folder = shapes_folders[i],
    mode = digitization_mode
  )
  
  all_jobs <- rbind(all_jobs, prep_result$jobs)
  all_temp_files <- c(all_temp_files, prep_result$temp_files)
  all_temp_dirs  <- c(all_temp_dirs, prep_result$temp_dir)
}

jobs <- all_jobs
temp_files <- all_temp_files
temp_dirs  <- unique(all_temp_dirs)

# Check if any jobs to process
if (nrow(jobs) == 0) {
  cat("No images to digitize in", digitization_mode, "mode.\n")
  if (digitization_mode == "skip") {
    cat("  - All images have already been digitized\n")
    cat("  - Use mode='review' to re-examine digitized images\n")
  } else {
    cat("  - No digitized images found to review\n")
    cat("  - Use mode='skip' to digitize new images\n")
  }
} else {
  # Randomize order
  jobs <- jobs[sample(nrow(jobs)), , drop = FALSE]
  
  regular_jobs <- 0
  rep_jobs <- 0
  imaging_rep_jobs <- 0
  
  for (i in 1:nrow(jobs)) {
    shape_file <- jobs$shape[i]
    if (grepl("_rep\\.txt$", shape_file)) {
      rep_jobs <- rep_jobs + 1
    } else if (grepl("2\\.txt$", shape_file)) {
      imaging_rep_jobs <- imaging_rep_jobs + 1
    } else if (grepl("1\\.txt$", shape_file)) {
      regular_jobs <- regular_jobs + 1
    }
  }
  
  cat("Launching StereoMorph for", wing_label, "of Bombus", paste(species, collapse = " + "), "\n")
  cat("Total images to digitize:", nrow(jobs), "\n")
  if (regular_jobs > 0) cat("  - Regular digitizations:", regular_jobs, "\n")
  if (rep_jobs > 0) cat("  - Digitization replicates:", rep_jobs, "\n")
  if (imaging_rep_jobs > 0) cat("  - Imaging replicates:", imaging_rep_jobs, "\n")
  cat("\n")
  
  # Launch StereoMorph
  digitizeImages(
    image.file    = jobs$image,
    shapes.file   = jobs$shape,
    landmarks.ref = lm_file
  )
  
  cat("Digitization complete!\n")
  
  # Clean up temporary files and directories
  if (length(temp_files) > 0) cleanup_temp_files(temp_files)
  
  for (td in temp_dirs) {
    if (dir.exists(td) && length(list.files(td)) == 0) {
      unlink(td, recursive = TRUE)
    }
  }
  
  # Summary
  cat("\nSummary:\n")
  cat("- Output folder:", shapes_folder, "\n")
}
