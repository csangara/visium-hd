library(tidyverse)
library(patchwork)
library(Seurat)

celltype_colors_df <- read.csv("data/ABA_metadata/WMB-taxonomy/20231215/views/cluster_to_cluster_annotation_membership_color.csv") %>% 
  inner_join(read.csv("data/ABA_metadata/WMB-taxonomy/20231215/views/cluster_to_cluster_annotation_membership_pivoted.csv"),
             by = "cluster_alias") %>% 
  mutate(class = str_replace_all(class, "[- ]", ""),
         class_color = ifelse(class_color == "#FFFB46","#DBD721", class_color)) %>% 
  distinct(class, class_color) %>% 
  arrange(class)
celltype_colors <-  celltype_colors_df$class_color %>% setNames(celltype_colors_df$class)

parcellations <- inner_join(read.csv("data/ABA_metadata/Allen-CCF-2020/20230630/views/parcellation_to_parcellation_term_membership_color.csv"),
                            read.csv("data/ABA_metadata/Allen-CCF-2020/20230630/views/parcellation_to_parcellation_term_membership_acronym.csv"),
                            by="parcellation_index")
division_colors_df <- parcellations %>% distinct(division, division_color) %>% 
  # If color is #FFFFFF or #CCCCCC, change division to "fiber tracts/VS"
  mutate(division = ifelse(division_color %in% c("#FFFFFF", "#CCCCCC"),
                           "fiber tracts/VS", division)) %>% 
  distinct(division, division_color) %>% 
  filter(division_color != "#FFFFFF") %>% 
  # Change all #AAAAAA colors to #CCCCCC
  mutate(division_color = ifelse(division_color == "#AAAAAA",
                                 "#CCCCCC", division_color))
division_colors <- setNames(division_colors_df$division_color, division_colors_df$division)
proper_celltype_names <- read.csv("data/ABA_metadata/WMB-taxonomy/20231215/views/cluster_to_cluster_annotation_membership_pivoted.csv") %>% 
  pull(class) %>% unique %>% sort %>% 
  setNames(names(celltype_colors))
proper_region_names <- c("Isocortex" = "Isocortex", "OLF" = "Olfactory areas", "HPF" = "Hippocampal formation", "CTXsp" = "Cortical subplate",
                         "STR" = "Striatum", "PAL" = "Pallidum", "TH" = "Thalamus", "HY" = "Hypothalamus")
sections_oi <- c("C57BL6J-638850.40","C57BL6J-638850.50")