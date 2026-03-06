library(Seurat)
library(tidyverse)

bin_size <- 32 # 8, 16, or 32

# Fill the text with 0
bin_size_str <- sprintf("%03dum", bin_size)

dataset <- "" #or "_caw009" or ""
data_path <- paste0("data/Visium_HD_Liver", toupper(dataset), "/")
plot_path <- paste0("visium_hd_liver", dataset, "/plots/")

visium_obj <- readRDS(paste0(data_path, "Visium_HD_Liver", toupper(dataset), "_",
                             bin_size_str, ".rds"))
dim(visium_obj) # 19059 genes x 462269 spots

p_ncount <- SpatialFeaturePlot(visium_obj, paste0("nCount_Spatial.", bin_size_str),
                               image.alpha=0, stroke=NA) +
  labs(fill = "nCount") +
  theme(legend.position = "right",
        legend.key.height = unit(3, "cm"),
        legend.key.width = unit(2, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))
#p_ncount
ggsave(paste0(plot_path, "spatial_nCount/spatialfeatureplot_nCount_", bin_size_str, ".png"),
       p_ncount,
       width = 40, height = 30, bg = "white") 

p_nfeatures <- SpatialFeaturePlot(visium_obj, paste0("nFeature_Spatial.", bin_size_str),
                               image.alpha=0, stroke=NA) +
  labs(fill = "nFeature") +
  theme(legend.position = "right",
        legend.key.height = unit(3, "cm"),
        legend.key.width = unit(2, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))
ggsave(paste0(plot_path, "spatial_nFeature/spatialfeatureplot_nFeature_", bin_size_str, ".png"),
       p_nfeatures,
       width = 40, height = 30, bg = "white")


# As these plots are a bit too huge, divide them into four quadrants
coords <- GetTissueCoordinates(visium_obj)

# Get quadrants
x_mid <- mean(c(max(coords$x), min(coords$x)))
y_mid <- mean(c(max(coords$y), min(coords$y)))

# Label each cell per quadrant
coords <- coords %>% 
  mutate(quadrant = case_when(x > x_mid & y > y_mid   ~ "Q1",
                              x <= x_mid & y > y_mid  ~ "Q2",
                              x <= x_mid & y <= y_mid ~ "Q3",
                              TRUE ~ "Q4"))

visium_obj$quadrant <- coords$quadrant

for (i in 1:4){
  p_ncount <- SpatialFeaturePlot(visium_obj %>% .[, .$quadrant == paste0("Q", i)],
                                 paste0("nCount_Spatial.", bin_size_str),
                                 stroke=NA, image.alpha=0, pt.size.factor = 5) +
    labs(fill = "nCount") +
    theme(legend.position = "right")
  
  #p_ncount
  ggsave(paste0(plot_path, "spatial_nCount/spatialfeatureplot_nCount_Q", i, "_", bin_size_str, ".png"),
         p_ncount,
         width = 10, height = 7, bg = "white") 
}

for (i in 1:4){
  p_nfeatures <- SpatialFeaturePlot(visium_obj %>% .[, .$quadrant == paste0("Q", i)],
                                    paste0("nFeature_Spatial.", bin_size_str),
                                    stroke=NA, image.alpha=0, pt.size.factor = 5) +
    labs(fill = "nFeature") +
    theme(legend.position = "right")
  
  ggsave(paste0(plot_path, "spatial_nFeature/spatialfeatureplot_nFeature_Q", i, "_", bin_size_str, ".png"),
         p_nfeatures,
         width = 10, height = 7, bg = "white")
}
