library(Seurat)
library(tidyverse)
library(patchwork)

## COLOR PALETTES ##
color_palette <- c("Hepatocytes" = "#B4B5B5FF",
                   "CentralVeinEndothelialcells" = "#FED8B1FF",
                   "LSECs" = "#FBB05FFF",
                   "PortalVeinEndothelialcells" = "#CC7722FF",
                   "LymphaticEndothelialcells" = "#8F4716FF",
                   "Cholangiocytes" = "#C61B84FF",
                   "HsPCs" = "#F19FC3FF",
                   "Stellatecells" = "#A31A2AFF",
                   "Mesothelialcells" = "#D0110BFF",
                   "Fibroblasts" = "#E45466FF",
                   "CapsularFibroblasts" = "#D46F6CFF",
                   "Kupffercells" = "#5DA6DBFF",
                   "MonocytesMonocytederivedcells" = "#a3daf3",
                   "cDC1s" = "#893A86FF",
                   "cDC2s" = "#893A86FF",
                   "pDCs" = "#893A86FF",
                   "MigcDCs" = "#893A86FF",
                   "Bcells" = "#9C7EBAFF",
                   "NKcells" = "#4A6E34FF",
                   "Tcells" = "#3AB04AFF",
                   "ILC1s" = "#A3D7BAFF",
                   "Basophils" = "#191919",
                   "Neutrophils" = "#727272")

color_palette_vizgen <- c("Hepatocytes" = "#B4B5B5FF",
                          "Endothelialcells" = "#dca754",
                          "Cholangiocytes" = "#C61B84FF",
                          "Stromalcells" = "#79151e",
                          "Kupffercells" = "#5DA6DBFF",
                          "Otherimmunecells" = "#893A86FF",
                          "Bcells" = "#9C7EBAFF",
                          "Unknown" = "#191919")

## FUNCTIONS ##
add_space_celltype <- function(celltype_strs) {
  celltype_strs %>%
    # Add space before the word "cells"
    stringr::str_replace_all("([A-Za-z]+)(cells)$", "\\1 \\2") %>%
    # Add space before the word "immune"
    stringr::str_replace_all("([A-Za-z]+)(immune)", "\\1 \\2")
}

celltype_labeller <- function(celltype_str){
  celltype_str %>% stringr::str_replace_all("Endothelialcells", "ECs") %>%
    stringr::str_replace_all("MonocytesMonocytederivedcells", "Mono & Mono-derived cells") %>%
    stringr::str_replace_all("([A-Za-z]+)(cells)$", "\\1 \\2") %>%
    stringr::str_replace_all("([a-z]{2,})([A-Z])", "\\1 \\2")
}

group_annot_to_vizgen <- function(data, input_col, output_col){
  mutate(data, !!output_col := case_when(grepl("Mesothelial|Stellate|VSMC|Fibroblasts", {{input_col}}) ~ "Stromalcells",
                                         grepl("Tcells|NKcells|Basophils|DC|Monocytes|ILC1s|Neutrophils|HsPCs", {{input_col}}) ~ "Otherimmunecells",
                                         grepl("Endothelialcells|ECs", {{input_col}}) ~ "Endothelialcells",
                                         # Others stay the same
                                         TRUE ~ {{input_col}})
  )
}


## THEMES ##
theme_histogram <- theme(strip.background = element_rect(linewidth = 0.5),
                         strip.text.y.right = element_text(size=7),
                         strip.text.x.top = element_text(hjust=0.5, size=7),
                         panel.spacing.y = unit(5, "mm"),
                         panel.grid.major.x = element_line(color="gray90", linewidth=0.1),
                         panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, color="black",linewidth = 0.25),
                         axis.title.x = element_blank(),
                         axis.title.y = element_text(size=6, margin = margin(r=7)),
                         axis.text.y = element_text(size=5),
                         axis.text.x = element_text(size=5),
                         axis.line.y.left = element_line(linewidth=0.25),
                         axis.line.x.bottom = element_line(linewidth=0.1),
                         axis.ticks = element_line(linewidth=0.1),
                         plot.title = element_text(size=7, face="bold"))

theme_barplot <-theme(panel.grid.minor.x = element_blank(),
                      panel.grid.major.x = element_blank(),
                      panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.25),
                      axis.ticks.x = element_blank(),
                      axis.ticks.y = element_line(linewidth=0.25),
                      axis.text.x = element_blank(),
                      axis.text.y = element_text(size=5),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(size=6),
                      plot.title = element_text(size=7, face="bold"),
                      legend.position = "bottom",
                      legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
                      legend.key.size = unit(0.4, "cm"),
                      legend.spacing.x = unit(0.4, "cm"),
                      legend.key.spacing.x = unit(0.3, "cm"),
                      legend.title = element_text(size=6, hjust=0.5, margin=margin(b=5)),
                      legend.text = element_text(size=6, margin = margin(l=3)))

theme_barplot_facet <- theme(panel.grid.major.y = element_blank(),
                             panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.5),
                             strip.background = element_blank(),
                             axis.ticks.y = element_blank(),
                             axis.ticks.x = element_line(linewidth=0.25),
                             axis.title.y = element_blank(),
                             axis.title.x = element_text(size=6),
                             plot.title = element_text(size=7, face="bold"),
                             legend.box.margin = margin(t = 0, r = 0, b = 0, l = 5),
                             legend.key.size = unit(0.5, "cm"),
                             legend.title = element_text(size=6),
                             legend.text = element_text(size=6),
                             legend.position = "bottom",
                             legend.direction = "vertical")

## OTHER VALUES ##
plot_path <- paste0("visium_hd_liver_combined/plots/")
celltype_order <- names(color_palette)
bin_sizes <- c(8, 16, 32)
bin_size <- 8
bin_size_strs <- sprintf("%03dum", bin_sizes)
