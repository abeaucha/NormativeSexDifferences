# label_human_expr.R
# Author: Antoine Beauchamp
# Edited: June 22nd, 2023
#
# 

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))


# Paths ----------------------------------------------------------------------

# Directories
expr_dir <- "data/human/expression"
atlas_dir <- "data/human/atlas"

# Microarray sample expression matrix
expr_file <- "HumanExpressionMatrix_samples_pipeline_abagen.csv"
expr_file <- file.path(expr_dir, expr_file)

# Microarray sample annotations
coordinates_file <- "AHBA_microarray_coordinates_mni.csv"
coordinates_file <- file.path(expr_dir, coordinates_file)

# Microarray sample annotations
annotations_file <- "AHBA_microarray_sample_annotations.csv"
annotations_file <- file.path(expr_dir, annotations_file)

# Microarray sample metadata
metadata_file <- "SampleInformation_pipeline_abagen.csv"
metadata_file <- file.path(expr_dir, metadata_file)

# Atlas definitions
defs_files <- c("glasser_defs.csv",
                "hypothalamus_defs.csv",
                "subcortical_defs.csv",
                "brainstem_defs.csv",
                "hippocampus_amygdala_defs.csv",
                "hippocampus_amygdala_defs.csv")
defs_files <- file.path(atlas_dir, defs_files)

# Atlas labels
labels_files <- c("glasser_labels.mnc",
                  "hypothalamus_labels.mnc",
                  "subcortical_labels.mnc",
                  "brainstem_labels.mnc",
                  "hippocampus_amygdala_labels_left.mnc",
                  "hippocampus_amygdala_labels_right.mnc")
labels_files <- file.path(atlas_dir, labels_files)


# Main -----------------------------------------------------------------------

# Import expression matrix
expr <- as_tibble(data.table::fread(expr_file, header = TRUE))

# Import microarray coordinates
coordinates <- read_csv(coordinates_file, show_col_types = FALSE)

# Import microarray sample annotations
annotations <- read_csv(annotations_file, show_col_types = FALSE)

#  Import microarray sample metadata
metadata <- read_csv(metadata_file, show_col_types = FALSE)

# Filter metadata for good samples
metadata <- metadata[annotations[["keep"]],]

# Export filtered metadata
outfile <- "human_sample_information.csv"
outfile <- file.path(expr_dir, outfile)
write_csv(x = metadata, file = outfile)

# Transpose expression matrix
expr <- expr %>%
  column_to_rownames("Gene") %>%
  as.matrix() %>%
  t() %>% 
  as_tibble()

# Filter microarray data for good samples
expr <- expr[annotations[["keep"]],]

# Initialize a tibble to store atlas regions
regions <- matrix(nrow = nrow(coordinates), 
                  ncol = length(defs_files))
colnames(regions) <- c("glasser", "hypothalamus",
                       "subcortical", "brainstem", 
                       "left_hipp_amyg", "right_hipp_amyg")
regions <- as_tibble(regions)

# Iterate over atlases and microarray samples
for (j in 1:ncol(regions)) {
  
  message(paste("Processing atlas:", colnames(regions)[j]))
  
  # Import atlas definitions
  defs <- read_csv(defs_files[j], show_col_types = FALSE)
  
  # Get atlas labels at microarray locations
  for (i in 1:nrow(coordinates)) {
    x <- coordinates[[i, "x"]]
    y <- coordinates[[i, "y"]]
    z <- coordinates[[i, "z"]]
    label <- mincGetWorldVoxel(filenames = labels_files[j], 
                               v1 = x, v2 = y, v3 = z)
    label <- round(label)
    if (label != 0) {
      regions[i,j] <- defs[label == defs[["label"]], "name"][[1]]
    }
  }
}
colnames(regions) <- str_c("Region_", colnames(regions))

# Label expression matrix with ROIs
expr <- bind_cols(expr, regions)

# Export labelled expression matrix
outfile <- "human_expression_matrix.csv"
outfile <- file.path(expr_dir, outfile)
data.table::fwrite(x = expr, file = outfile)
