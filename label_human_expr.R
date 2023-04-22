# label_human_expr.R
# Author: Antoine Beauchamp
# Edited: April 22nd, 2023
#
# 

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))


# Directories ----------------------------------------------------------------

expr_dir <- "data/human/expression"
atlas_dir <- "data/human/atlas"


# Main -----------------------------------------------------------------------

#Import expression matrix
expr_file <- "HumanExpressionMatrix_samples_pipeline_abagen.csv"
expr_file <- file.path(expr_dir, expr_file)
expr <- as_tibble(data.table::fread(expr_file, header = TRUE))

#Import microarray coordinates
coordinates_file <- "AHBA_microarray_coordinates_mni.csv"
coordinates_file <- file.path(expr_dir, coordinates_file)
coordinates <- read_csv(coordinates_file, show_col_types = FALSE)

#Import microarray sample annotations
annotations_file <- "AHBA_microarray_sample_annotations.csv"
annotations_file <- file.path(expr_dir, annotations_file)
annotations <- read_csv(annotations_file, show_col_types = FALSE)

#Import atlas definitions
defs_file <- "glasser_hypothalamus_defs.csv"
defs_file <- file.path(atlas_dir, defs_file)
defs <- read_csv(defs_file, show_col_types = FALSE)

#Labels file
labels_file <- "glasser_hypothalamus_labels.mnc"
labels_file <- file.path(atlas_dir, labels_file)

#Transpose expression matrix
expr <- expr %>%
  column_to_rownames("Gene") %>%
  as.matrix() %>%
  t() %>% 
  as_tibble()

#Filter microarray data for good samples
expr <- expr[annotations[["keep"]],]

#Get atlas labels at microarray locations
structs <- character(nrow(coordinates))
for (i in 1:nrow(coordinates)) {
  x <- coordinates[[i, "x"]]
  y <- coordinates[[i, "y"]]
  z <- coordinates[[i, "z"]]
  label <- mincGetWorldVoxel(filenames = labels_file, v1 = x, v2 = y, v3 = z)
  label <- round(label)
  if (label != 0) {
    structs[i] <- defs[label == defs[["label"]], "name"][[1]]
  }
}

#Label expression matrix with ROIs
expr[["Region"]] <- structs

#Export labelled expression matrix
outfile <- "human_expression_matrix.csv"
outfile <- file.path(expr_dir, outfile)
data.table::fwrite(x = expr, file = outfile)