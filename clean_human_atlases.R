# clean_human_atlases.R
# Author: Antoine Beauchamp
# Edited: May 15th, 2023


# Packages -------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))


# Paths ----------------------------------------------------------------------

# Directories 
atlas_dir <- "data/human/atlas/"

# Glasser atlas definitions
defs_glasser <- "HCPex.csv"
defs_glasser <- file.path(atlas_dir, defs_glasser)

# Hypothalamus atlas definitions
defs_hypo <- "Volumes_names-labels.csv"
defs_hypo <- file.path(atlas_dir, defs_hypo)

# Subcortical atlas definitions
defs_subcx <- "FreeSurferColorLUT.csv"
defs_subcx <- file.path(atlas_dir, defs_subcx)

# Brainstem atlas definitions
defs_bs <- "compressionLookupTable_brainstem.csv"
defs_bs <- file.path(atlas_dir, defs_bs)

# Hippocampus/amygdala atlas definitions
defs_hippamyg <- "compressionLookupTable_HP_Amyg.csv"
defs_hippamyg <- file.path(atlas_dir, defs_hippamyg)


# Main -----------------------------------------------------------------------

# Clean Glasser atlas definitions
outfile <- "glasser_defs.csv"
outfile <- file.path(atlas_dir, outfile)
read_csv(defs_glasser, show_col_types = FALSE) %>% 
  select(hemisphere = Hemisphere,
         name = Label_name,
         label = Label_number...1) %>%
  mutate(hemisphere = ifelse(hemisphere == "L", "left", "right")) %>% 
  unite(col = "name", hemisphere, name, sep = " ") %>% 
  write_csv(file = outfile)

# Clean hypothalamus atlas definitions
outfile <- "hypothalamus_defs.csv"
outfile <- file.path(atlas_dir, outfile)
read_csv(defs_hypo, show_col_types = FALSE) %>% 
  select(hemisphere = Hemisphere,
         name = Name,
         label = Label) %>% 
  unite(col = "name", hemisphere, name, sep = " ") %>% 
  write_csv(file = outfile)

# Clean subcortical atlas definitions
outfile <- "subcortical_defs.csv"
outfile <- file.path(atlas_dir, outfile)
read_csv(defs_subcx, show_col_types = FALSE) %>% 
  select(name, label) %>% 
  write_csv(file = outfile)

# Clean brainstem atlas definitions
outfile <- "brainstem_defs.csv"
outfile <- file.path(atlas_dir, outfile)
read_csv(defs_bs, show_col_types = FALSE) %>% 
  select(name = label_name,
         label = Not_sure) %>% 
  write_csv(outfile)

# Clean hippocampus/amygdala atlas definitions
outfile <- "hippocampus_amygdala_defs.csv"
outfile <- file.path(atlas_dir, outfile)
read_csv(defs_hippamyg, show_col_types = FALSE) %>% 
  select(name = label_name,
         label = Not_sure) %>% 
  mutate(name = str_remove(name, "^Left-")) %>% 
  write_csv(file = outfile)
