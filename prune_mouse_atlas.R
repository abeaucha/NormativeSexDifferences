# prune_mouse_atlas.R
# Author: Antoine Beauchamp
# Edited: April 22nd, 2023
#
#
#

#Packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.tree))
suppressPackageStartupMessages(library(RMINC))

#Functions
source("src/tree_tools.R")

#Leaf nodes for pruned lateral tree
leaf_nodes_lateral <- c("Piriform cortex",
                        "Subiculum",
                        "Entorhinal area",
                        "Field CA1",
                        "Field CA3",
                        "Dentate gyrus, molecular layer",
                        "Medial amygdalar nucleus",
                        "Cortical subplate",
                        "Anterior cingulate area",
                        "Retrosplenial area",
                        "Primary auditory area",
                        "Agranular insular area",
                        "Perirhinal area",
                        "Primary motor area",
                        "Orbital area, ventrolateral part",
                        "Primary somatosensory area",
                        "Primary visual area",
                        "Temporal association areas",
                        "Posterior parietal association areas",
                        "Pallidum",
                        "Bed nuclei of the stria terminalis",
                        "Striatum ventral region",
                        "Caudoputamen",
                        "Thalamus",
                        "Mammillary body",
                        "Medial preoptic nucleus",
                        "Cerebellar cortex")

#Leaf nodes for pruned midline tree
leaf_nodes_midline <- c("Midbrain", "Medulla", "Pons")

#Load mouse neuroanatomical tree
treefile <- "data/mouse/expression/MouseExpressionTree_DSURQE.RData"
load(treefile)
tree <- Clone(treeMouseExpr)
rm(treeMouseExpr)

#Remove vermal regions, white matter and ventricles
tree <- prune_tree(tree = tree, 
                   nodes = c("Vermal regions", 
                             "fiber tracts",
                             "ventricular systems"), 
                   method = "at", inplace = FALSE)

#Create left hemispheric tree
tree_left <- prune_tree(tree = tree, 
                        nodes = leaf_nodes_midline,
                        method = "at")
tree_left <- lateralize_tree(tree = tree_left, hemisphere = "left")
tree_left <- prune_tree(tree = tree_left,
                        nodes = paste("left", leaf_nodes_lateral),
                        method = "below",
                        inplace = FALSE)

#Create right hemispheric tree
tree_right <- prune_tree(tree = tree, 
                        nodes = leaf_nodes_midline,
                        method = "at")
tree_right <- lateralize_tree(tree = tree_right, hemisphere = "right")
tree_right <- prune_tree(tree = tree_right,
                        nodes = paste("right", leaf_nodes_lateral),
                        method = "below",
                        inplace = FALSE)

#Create midline tree
tree_midline <- prune_tree(tree = tree,
                           nodes = leaf_nodes_midline,
                           method = "below",
                           remove = TRUE)

#Import DSURQE labels
atlas_dir <- "data/mouse/atlas"
labels_file <- "DSURQE_CCFv3_labels_200um.mnc"
labels_file <- file.path(atlas_dir, labels_file)
labels <- mincGetVolume(labels_file)

#Create left, right and midline aggregated label images
labels_left <- hanatToAtlas(tree_left, labelVolume = mincArray(labels))
labels_right <- hanatToAtlas(tree_right, labelVolume = mincArray(labels))
labels_midline <- hanatToAtlas(tree_midline, labelVolume = mincArray(labels))

#Empty atlas labels image
labels_pruned <- numeric(length(labels)) 

#Add left labels
ind_left <- labels_left > 0
labels_pruned[ind_left] <- labels_left[ind_left]

#Add right labels
ind_right <- labels_right > 0
labels_pruned[ind_right] <- labels_right[ind_right]

#Add midline labels
ind_midline <- labels_midline > 0
labels_pruned[ind_midline] <- labels_midline[ind_midline]

#Copy attributes from DSURQE
attributes(labels_pruned) <- attributes(labels)

#Export pruned atlas labels
labels_pruned_file <- str_replace(basename(labels_file), ".mnc", "_pruned.mnc")
labels_pruned_file <- file.path(atlas_dir, labels_pruned_file)
mincWriteVolume(buffer = labels_pruned,
                output.filename = labels_pruned_file,
                like.filename = labels_file,
                clobber = TRUE)

#Get left atlas definitions
defs_left <- hanatToAtlasDefs(tree_left) %>% 
  rename(name = Structure,
         label = Label)

#Get right atlas definitions
defs_right <- hanatToAtlasDefs(tree_right) %>% 
  rename(name = Structure,
         label = Label)

#Get midline atlas definitions
defs_midline <- hanatToAtlasDefs(tree_midline) %>% 
  rename(name = Structure,
         label = Label)

#Combine atlas defs
defs_pruned <- bind_rows(defs_left, defs_right, defs_midline)

#Export definitions
defs_pruned_file <- "DSURQE_defs_pruned.csv"
defs_pruned_file <- file.path(atlas_dir, defs_pruned_file)
write_csv(x = defs_pruned, file = defs_pruned_file)

