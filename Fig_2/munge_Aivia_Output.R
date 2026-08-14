library(dplyr)
library(readr)

##### Config

# Set the main directory containing all image folders
main_dir <- "C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Microscopy Analysis/20260129a_FigureRNAtreat/2026-01-29-16-30-54 IF_analysis/"   # UPDATE THIS PATH

# Optional: specify a particular measurement file name
# If NULL, will use the first CSV file found in each Measurements folder
measurement_filename <- "Measurements_Set_Meshes.csv"  # or set to "your_file.csv"

###### Process folders

# Get all subdirectories in the main folder
folders <- list.dirs(main_dir, recursive = FALSE, full.names = TRUE)

cat("Found", length(folders), "folders to process\n\n")

# Initialize empty list to store data frames
results_list <- list()

# Loop through each folder
for (folder in folders) {
  folder_name <- basename(folder)
  cat("Processing folder:", folder_name, "\n")
  
  # Find .aivia.tif file to get image name
  aivia_files <- list.files(folder, pattern = "\\.aivia\\.tif$", full.names = FALSE)
  
  if (length(aivia_files) == 0) {
    cat("  Skipping - no .aivia.tif file found\n\n")
    next
  }
  
  if (length(aivia_files) > 1) {
    warning(paste("Multiple .aivia.tif files found in:", folder_name, 
                  "- using first one:", aivia_files[1]))
  }
  
  # Extract image name (remove .aivia.tif extensions)
  # Remove everything after the last two periods
  image_name <- sub("\\.[^.]+\\.[^.]+$", "", aivia_files[1])
  cat("  Image name:", image_name, "\n")
  
  # Check if Measurements folder exists
  measurements_dir <- file.path(folder, "Measurements")
  
  if (!dir.exists(measurements_dir)) {
    cat("  Skipping - Measurements folder not found\n\n")
    next
  }
  
  # Determine which CSV file to read
  if (is.null(measurement_filename)) {
    # Find CSV files in Measurements folder
    csv_files <- list.files(measurements_dir, pattern = "\\.csv$", full.names = TRUE)
    
    if (length(csv_files) == 0) {
      cat("  Skipping - no CSV files found in Measurements folder\n\n")
      next
    }
    
    if (length(csv_files) > 1) {
      warning(paste("Multiple CSV files found in:", measurements_dir,
                    "- using first one:", basename(csv_files[1])))
    }
    
    measurement_path <- csv_files[1]
  } else {
    # Use specified filename
    measurement_path <- file.path(measurements_dir, measurement_filename)
    
    if (!file.exists(measurement_path)) {
      cat("  Skipping - specified measurement file not found\n\n")
      next
    }
  }
  
  # Read the measurement table
  tryCatch({
    df <- read_csv(measurement_path, show_col_types = FALSE)
    cat("  Read", nrow(df), "rows\n")
    
    # Add image identifier column as the first column
    df <- df %>%
      mutate(Image = image_name, .before = 1)
    
    # Add to results list
    results_list[[folder_name]] <- df
    cat("  Successfully processed\n\n")
    
  }, error = function(e) {
    cat("  Error reading file:", e$message, "\n\n")
  })
}

###### Combine results 

# Check if any data was collected
if (length(results_list) == 0) {
  stop("No measurement tables were successfully read")
}

# Combine all tables with rbind
master_table <- bind_rows(results_list)




###### Save
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Microscopy Analysis/20260129a_FigureRNAtreat/2026-01-29-16-30-54 IF_analysis/")
write_csv(master_table, "combined_measurements_2026-01-29-16-30-54 IF_analysis_A.csv")
