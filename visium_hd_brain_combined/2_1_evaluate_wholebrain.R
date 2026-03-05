library(tidyverse)
library(Seurat)
source("visium_hd_brain_combined/0_brain_init.R")
source("scripts/scatterbarplot_function.R")

# Read RCTD files
deconv_props_list <- list()
visium_objs_list <- list()
for (ds in c("ffpe", "fresh_frozen")){
  ds_ext <- ifelse(ds == "ffpe", "", "_FF")
  
  visium_obj <- readRDS(paste0("data/Visium_HD_MouseBrain", ds_ext, "/Visium_HD_MouseBrain", ds_ext, "_008um.rds"))
  deconv_props <- read.table(paste0("visium_hd_brain", tolower(ds_ext), "/Visium_HD_MouseBrain", ds_ext,
                                    "_008um/proportions_rctd_Visium_HD_MouseBrain", ds_ext, "_008um"),
                             header = TRUE)
  
  removed_spots <- Cells(visium_obj)[which(visium_obj@meta.data[,"nCount_Spatial.008um"] < 100)]
  
  stopifnot(length(removed_spots) + dim(deconv_props)[1] == dim(visium_obj)[2])
  visium_obj_subset <- visium_obj[, !(colnames(visium_obj) %in% removed_spots)]
  rownames(deconv_props) <- colnames(visium_obj_subset)
  
  visium_objs_list[[ds]] <- visium_obj_subset
  deconv_props_list[[ds]] <- deconv_props
  
  rm(visium_obj, visium_obj_subset, deconv_props, removed_spots)
}

deconv_props_df <- lapply(1:2, function(i){
  deconv_props_list[[i]] %>% tibble::rownames_to_column("spot") %>% 
    pivot_longer(-spot, names_to = "celltype", values_to = "proportion") %>% 
    mutate(celltype = str_replace_all(celltype, "^X", ""),
           source = names(deconv_props_list)[i])
}) %>% bind_rows()
  

first_run <- FALSE
ncelltypes <- 1
deconv_props_df_shapes <- lapply(unique(deconv_props_df$source), function(src) {
  
  square_size <- visium_objs_list[[src]]@images[["slice1.008um"]]@scale.factors$spot
  
  if (first_run){
    deconv_props_roi <- deconv_props_df %>%
      filter(source == src) %>%
      filter(proportion > 0) %>% 
      # merge with coordinates
      left_join(GetTissueCoordinates(visium_objs_list[[src]]),
                by = c("spot" = "cell"))
    
    if (src == "fresh_frozen"){
      deconv_props_roi <- deconv_props_roi %>% rename(x=y, y=x)
    }

    deconv_props_df_shapes <- plot_scatterbar(deconv_props_roi,
                                              square_size,
                                              ncelltypes=ncelltypes,
                                              return_df = TRUE)$df
    
    deconv_props_df_shapes <- deconv_props_df_shapes %>% 
      # Remove X from celltype 
      mutate(celltype = gsub("^X", "", celltype)) %>% 
      mutate(celltype = factor(celltype, levels = names(celltype_colors)))
    
    saveRDS(deconv_props_df_shapes, paste0("visium_hd_brain_combined/rds/deconv_props_shapes_",
                                             tolower(src), "_ncelltypes", ncelltypes, "_008um.rds"))
    
  } else {
    deconv_props_df_shapes <- readRDS(paste0("visium_hd_brain_combined/rds/deconv_props_shapes_",
                                             tolower(src), "_ncelltypes", ncelltypes, "_008um.rds"))
  }
  deconv_props_df_shapes$source <- src
  deconv_props_df_shapes
  
}) %>% bind_rows()


p_deconv_celltypes <- lapply(c("ffpe", "fresh_frozen"), function(src){
  ggplot(deconv_props_df_shapes %>% filter(source == src),
         aes(x = x, y = y)) +
    geom_polygon(aes(fill = celltype, group = group), linewidth = 0, show.legend = FALSE) +
    scale_fill_manual(values = celltype_colors) +
    {if (src == "fresh_frozen") scale_x_reverse() } +
    theme_void(base_size=7) +
    scale_y_reverse() +
    coord_fixed(ratio = 1)
})


# Read merfish
merfish_cells <- readRDS("visium_hd_brain_combined/rds/merfish_cells_3datasets_combined.rds")
merfish_cells %>% head
merfish_cells$brain_section_label %>% table

merfish_cells_section <- merfish_cells %>% 
  filter(brain_section_label %in% sections_oi[1], x < 5.5,
         !division %in% c("brain-unassigned", "unassigned")) %>% 
  mutate(class = str_replace_all(class, "[- ]", ""))

p_merfish_celltype <- ggplot(merfish_cells_section %>%
                               filter(brain_section_label == sections_oi[1]),
       aes(y=y, x=x, color=class)) +
  geom_point(shape=16, size=0.4, stroke=0, show.legend=FALSE) +
  scale_y_reverse() + coord_fixed(ratio=1) +
  scale_color_manual(values = celltype_colors) +
  theme_void(base_size=7) +
  theme(legend.position = "bottom",
        axis.title = element_blank())


p_celltypes_combined <- p_merfish_celltype + p_deconv_celltypes[[1]] + p_deconv_celltypes[[2]]
p_celltypes_combined
ggsave("visium_hd_brain_combined/plots/merfish_vs_visiumhd_celltype_annotations.pdf",
         plot = p_celltypes_combined,
         width = 7, height = 5)
