# label_mouse_expr.R
# Author: Antoine Beauchamp
# Edited: April 22nd, 2023
#
# 

# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))


# Functions ------------------------------------------------------------------

#' Calculate a robust sigmoid transformation
#'
#' @param x (numeric vector) Input data to transform
#' @return (numeric vector) Transformed data
robust_sigmoid <- function(x) {
  IQR_x <- IQR(x)/(qnorm(0.75) - qnorm(0.25))
  y <- 1/(1 + exp(-1*(x - median(x))/IQR_x))
  return(y)
}


#' Calculate a min-max transformation
#'
#' @param x (numeric vector) Input data to transform
#' @return (numeric vector) Transformed data
min_max <- function(x){(x - min(x))/(max(x) - min(x))}


#' Calculate a scaled robust sigmoid transformation
#'
#' @param x (numeric vector) Input data to transform
#' @return (numeric vector) Transformed data
scaled_robust_sigmoid <- function(x) {min_max(robust_sigmoid(x = x))}


# Paths ----------------------------------------------------------------------

#Directories
expr_dir <- "data/mouse/expression/"
atlas_dir <- "data/mouse/atlas"

#Voxelwise expression matrix
expr_file <- "MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed.csv"
expr_file <- file.path(expr_dir, expr_file)

#Atlas definitions
defs_file <- "DSURQE_defs_pruned.csv"
defs_file <- file.path(atlas_dir, defs_file)

#Atlas labels
labels_file <- "DSURQE_CCFv3_labels_200um_pruned.mnc"
labels_file <- file.path(atlas_dir, labels_file)

#Mask
mask_file <- "coronal_200um_coverage_bin0.8.mnc"
mask_file <- file.path(atlas_dir, mask_file)


# Main -----------------------------------------------------------------------

#Import expression matrix
expr <- as_tibble(data.table::fread(expr_file, header = TRUE))

#Import atlas definitions
defs <- read_csv(defs_file, show_col_types = FALSE)

#Import atlas labels
labels <- round(mincGetVolume(labels_file))

#Import brain mask
mask <- mincGetVolume(mask_file)

#Normalize expression matrix
#Note: There's an implicitly transpose resulting from the second apply()
expr <- expr %>% 
  column_to_rownames("Gene") %>% 
  as.matrix() %>% 
  apply(MARGIN = 2, FUN = scaled_robust_sigmoid) %>% 
  apply(MARGIN = 1, FUN = scaled_robust_sigmoid) %>% 
  as_tibble()

#Apply mask to labels
labels_masked <- labels[mask > 0]

#Match voxel labels to definition
ind_match <- match(labels_masked, defs[["label"]])
structs <- defs[["name"]][ind_match]

#Label expression matrix with ROIs
expr[["Region"]] <- structs

#Export labelled expression matrix
outfile <- "mouse_expression_matrix.csv"
outfile <- file.path(expr_dir, outfile)
data.table::fwrite(x = expr, file = outfile)
