# stitch_human_atlases.R
# Author: Antoine Beauchamp
# Edited: April 22nd, 2023


# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RMINC))


# Directories ----------------------------------------------------------------

atlas_dir <- "data/human/atlas/"


# Paths ----------------------------------------------------------------------

#Glasser atlas definitions
defs_glasser <- "HCPex.csv"
defs_glasser <- file.path(atlas_dir, defs_glasser)

#Glasser atlas labels
labels_glasser <- "Glasser_atlas_in_9cSym.mnc"
labels_glasser <- file.path(atlas_dir, labels_glasser)

#Hypothalamus atlas definitions
defs_hypo <- "Volumes_names-labels.csv"
defs_hypo <- file.path(atlas_dir, defs_hypo)

#Hypothalamus atlas labels
labels_hypo <- "Hypothal_atlas_in_imcb9c.mnc"
labels_hypo <- file.path(atlas_dir, labels_hypo)

#MNI template
template_file <- "mni_icbm152_t1_tal_nlin_sym_09c.mnc"
template_file <- file.path(atlas_dir, template_file)


# Main -----------------------------------------------------------------------

#Import Glasser atlas and fix columns
defs_glasser <- read_csv(defs_glasser, show_col_types = FALSE)
defs_glasser <- defs_glasser %>% 
  select(hemisphere = Hemisphere,
         name = Label_name,
         label = Label_number...1) %>%
  mutate(hemisphere = ifelse(hemisphere == "L", "left", "right")) %>% 
  unite(col = "name", hemisphere, name, sep = " ") %>% 
  mutate(atlas = "glasser")

#Import hypothalamus atlas and fix columns
defs_hypo <- read_csv(defs_hypo, show_col_types = FALSE)
defs_hypo <- defs_hypo %>% 
  select(hemisphere = Hemisphere,
         name = Name,
         label = Label) %>% 
  unite(col = "name", hemisphere, name, sep = " ") %>% 
  mutate(atlas = "hypothalamus")

#Stitch atlas definitions 
defs_stitched <- bind_rows(defs_glasser, defs_hypo)
defs_stitched <- defs_stitched %>% 
  mutate(label_new = 1:nrow(.))

#Import MNI template
template <- mincGetVolume(template_file)

#Stitch atlases
atlases <- c("glasser", "hypothalamus")
label_files <- c(labels_glasser, labels_hypo)
labels_stitched <- numeric(length(template))
for (i in 1:length(atlases)) {
  labels <- mincGetVolume(label_files[i])
  defs <- defs_stitched %>% 
    filter(atlas == atlases[i])
  
  ind_match <- match(labels, defs[["label"]])
  labels_new <- defs[["label_new"]][ind_match]
  labels_new[is.na(labels_new)] <- 0
  
  ind_nonzero <- labels_new != 0
  labels_stitched[ind_nonzero] <- labels_new[ind_nonzero]
}
attributes(labels_stitched) <- attributes(labels)

#Export stitched atlas
labels_outfile <- "glasser_hypothalamus_labels.mnc"
labels_outfile <- file.path(atlas_dir, labels_outfile)
mincWriteVolume(buffer = labels_stitched, 
                output.filename = labels_outfile,
                like.filename = template_file,
                clobber = TRUE)

#Consolidate definitions
defs_stitched <- defs_stitched %>% 
  select(-label,
         -atlas, 
         label = label_new)

#Export definitions
defs_outfile <- "glasser_hypothalamus_defs.csv"
defs_outfile <- file.path(atlas_dir, defs_outfile)
write_csv(x = defs_stitched, file = defs_outfile)
