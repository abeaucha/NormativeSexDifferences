library(tidyverse)
library(RMINC)

atlas_dir <- "data/human/atlas/"

#Import Glasser atlas and fix columns
defs_glasser <- "HCPex.csv"
defs_glasser <- file.path(atlas_dir, defs_glasser)
defs_glasser <- read_csv(defs_glasser)
defs_glasser <- defs_glasser %>% 
  select(hemisphere = Hemisphere,
         name = Label_name,
         label = Label_number...1) %>%
  mutate(hemisphere = ifelse(hemisphere == "L", "left", "right")) %>% 
  unite(col = "name", hemisphere, name, sep = " ") %>% 
  mutate(atlas = "glasser")

#Import hypothalamus atlas and fix columns
defs_hypo <- "Volumes_names-labels.csv"
defs_hypo <- file.path(atlas_dir, defs_hypo)
defs_hypo <- read_csv(defs_hypo)
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
template_file <- "mni_icbm152_t1_tal_nlin_sym_09c.mnc"
template_file <- file.path(atlas_dir, template_file)
template <- mincGetVolume(template_file)

#Atlas label files
label_files <- c("Glasser_atlas_in_9cSym.mnc", "Hypothal_atlas_in_imcb9c.mnc")
label_files <- file.path(atlas_dir, label_files)

#Stitch atlases
atlases <- c("glasser", "hypothalamus")
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

