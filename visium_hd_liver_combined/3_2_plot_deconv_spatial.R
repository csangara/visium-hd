source("visium_hd_liver_combined/0_utils.R")
source("scripts/scatterbarplot_function.R")

first_run <- FALSE
rois <- list("caw009" = c(6600, 7100, 7200, 8200), # roi1
             "caw009" = c(9400, 9900, 11200, 12200), # roi2 (spp1)
             "sca002" = c(14000,15000, 12500, 14500), # roi1
             "sca002" = c(17500, 18500, 29000, 31000) # roi2 (spp1)
)


##### FIRST RUN: SAVE ROI OBJECTS #####

if (first_run){
for (dataset_name in c("sca002", "caw009")){
  cat("Processing dataset:", dataset_name, "\n")
  dataset <- ifelse(dataset_name == "caw009", "_caw009", "")
  data_path <- paste0("data/Visium_HD_Liver", toupper(dataset), "/")
  
  for (bin_size in bin_sizes){
    gc()
    bin_size_str <- sprintf("%03dum", bin_size)
    cat("Processing bin size:", bin_size_str, "\n")
    
    visium_obj <- readRDS(paste0(data_path, "Visium_HD_Liver", toupper(dataset), "_",
                                 bin_size_str, ".rds"))
    rois_oi <- rois[names(rois) == dataset_name]
    
    for (roi in rois_oi) {
      
      roi_name <- paste0("_", paste0(sprintf("%.1f", roi/1000), collapse="-"))
      cells_roi <- GetTissueCoordinates(visium_obj) %>%
        filter(x > roi[1] & x < roi[2] & y > roi[3] & y < roi[4]) %>% pull(cell)
      
      visium_obj_roi <- visium_obj %>% .[, colnames(.) %in% cells_roi]
      
      saveRDS(visium_obj_roi, paste0("visium_hd_liver_combined/rds/", "liver", toupper(dataset), "_",
                                     bin_size_str, "_ROI", roi_name, ".rds"))
    }
  }
}
}

##### FIRST RUN: SAVE HIRES IMAGE OBJECTS #####
# This draws the red rectangle over the ROI (Figure 4.5a, left)
if (first_run){
bin_size <- 16
bin_size_str <- sprintf("%03dum", bin_size)

for (dataset_name in c("sca002", "caw009")){
  dataset <- ifelse(dataset_name == "caw009", "_caw009", "")
  data_path <- paste0("data/Visium_HD_Liver", toupper(dataset), "/")
  rois_oi <- rois[names(rois) == dataset_name] %>% do.call(rbind, .) %>% 
    `colnames<-`(c("ymin", "ymax", "xmin", "xmax")) %>% 
    data.frame

  image_dir <- file.path(data_path, paste0("binned_outputs/square_", bin_size_str, "/spatial"))
  
  image_hires <- Read10X_Image(
    image_dir,
    image.name = "tissue_hires_image.png",
    assay = paste0("Spatial.", bin_size_str),
  )
  
  visium_obj_hires <- Load10X_Spatial(
    data_path,
    bin.size = bin_size,
    image = image_hires,
  )
  
  # For rectangle border: scale and shift y-axis so it starts from the top
  y_range <- range(GetTissueCoordinates(visium_obj_hires)$x)
  hires_scalefactor <- visium_obj_hires@images[[paste0("slice1.", bin_size_str)]]@scale.factors$hires
  
  rois_oi_scaled <- rois_oi %>% 
    mutate(xmin = xmin * hires_scalefactor,
           xmax = xmax * hires_scalefactor,
           ymin = (y_range[2] - ymin + y_range[1])*hires_scalefactor,
           ymax = (y_range[2] - ymax + y_range[1])*hires_scalefactor,
           label = paste0("ROI ", row_number()))

  # Compare this with the plot we get
  # cells_roi1 <- GetTissueCoordinates(visium_obj_hires) %>%
  #   filter(x > rois_oi[1,1] & x < rois_oi[1,2] & y > rois_oi[1,3] & y < rois_oi[1,4]) %>% pull(cell)
  # cells_roi2 <- GetTissueCoordinates(visium_obj_hires) %>%
  #   filter(x > rois_oi[2,1] & x < rois_oi[2,2] & y > rois_oi[2,3] & y < rois_oi[2,4]) %>% pull(cell)
  # cells_roi <- c(cells_roi1, cells_roi2)
  # SpatialDimPlot(visium_obj_hires,
  #                cells.highlight = cells_roi,
  #                image.scale="hires")

  p_roi <- SpatialDimPlot(visium_obj_hires, pt.size.factor = 0,
                          image.scale="hires") +
    geom_rect(data = rois_oi_scaled,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax),
              fill = NA, color = "red", linewidth = 0.5,
              inherit.aes = FALSE) +
    geom_text(data = rois_oi_scaled,
              aes(x = (xmin + xmax) / 2, y = ymin,
                  label = label),
              color = "red", size = 3, vjust = -1,
              inherit.aes = FALSE)
  
  ggsave(paste0(plot_path, "spatialdimplot", toupper(dataset),
                "_rois_hires.png"), p_roi,
         width = 8, height = 8, bg = "white")
  
  saveRDS(p_roi, paste0("visium_hd_liver_combined/rds/",
                        "spatialdimplot_ROIs", toupper(dataset), "_",
                 bin_size_str, ".rds"))
}
}

##### Figure 4.5 ####
dataset_name <- "caw009"
bin_size <- 16
bin_size_str <- sprintf("%03dum", bin_size)
features <- c("Glul", "Hal", "Spp1", "Clec4f")

for (dataset_name in c("sca002", "caw009")){

  dataset <- ifelse(dataset_name == "caw009", "_caw009", "")
  rois_oi <- rois[names(rois) == dataset_name]
  roi_names <- sapply(rois_oi, function(roi) paste0("_", paste0(sprintf("%.1f", roi/1000), collapse="-")))
  
  # Load the ROI dimplot
  p_roi <- readRDS(paste0("visium_hd_liver_combined/rds/",
                          "spatialdimplot_ROIs", toupper(dataset), "_",
                          bin_size_str, ".rds"))
  
  square_size <- readRDS(paste0("visium_hd_liver_combined/rds/","liver", toupper(dataset), "_",
                                bin_size_str, "_ROI", roi_names[1], ".rds"))@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot
  
  # Load the Visium ROI objects
  features_df <- lapply(1:length(roi_names), function(i) {
    roi_name <- roi_names[i]
    visium_obj_roi <- readRDS(paste0("visium_hd_liver_combined/rds/",
                                     "liver", toupper(dataset), "_",
                                     bin_size_str, "_ROI", roi_name, ".rds"))
    
    visium_obj_roi <- NormalizeData(visium_obj_roi)
    
    # Check if the object is loaded correctly
    print(paste("Loaded ROI object for", roi_name))
    
    features_df <- GetAssayData(visium_obj_roi, layer = "data",
                                assay = paste0("Spatial.", bin_size_str))[features,] %>% 
      t() %>% data.frame() %>% 
      rownames_to_column("cell") %>%
      inner_join(GetTissueCoordinates(visium_obj_roi), by = "cell") %>% 
      pivot_longer(cols = -c(cell, x, y), names_to = "feature", values_to = "value") %>% 
      mutate(roi = i)
    
    features_df
  }) %>% bind_rows()
  
  feature_roi <- paste0(rep(features, each=length(roi_names)), "_", 1:length(roi_names))
  
  # Feature plots
  p_features_list <- lapply (feature_roi, function(txt) {
    feature_name <- strsplit(txt, "_")[[1]][1]
    roi_i <- as.numeric(strsplit(txt, "_")[[1]][2])
  
    p <- ggplot(features_df %>% filter(feature == feature_name, roi == roi_i),
                aes(x = y, y = x)) +
      geom_tile(aes(fill = value), color="gray90", linewidth=0.2,
                height=square_size, width=square_size) +
      scale_fill_viridis_c(breaks=scales::breaks_extended(n=5),
                           limits = range(features_df$value),
                           name = "log1p counts") +
      labs(title = feature_name) +
      scale_y_reverse() +
      coord_fixed() +
      theme_void()
    
    p
  }) %>%  setNames(feature_roi)
  
  p_roi <- p_roi + theme(legend.position="none")
  
  binsize_roi <- paste0(rep(bin_sizes, each=length(roi_names)), "_", 1:length(roi_names))
  
  # Create the ggplot
  p_props_list <- lapply(binsize_roi, function(txt) {
    bs <- as.numeric(strsplit(txt, "_")[[1]][1])
    roi_i <- as.numeric(strsplit(txt, "_")[[1]][2])
    bin_size_str <- sprintf("%03dum", bs)
    
    visium_obj_roi <- readRDS(paste0("visium_hd_liver_combined/rds/",
                                     "liver", toupper(dataset), "_",
                                     bin_size_str, "_ROI", roi_names[roi_i], ".rds"))
    square_size <- visium_obj_roi@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot
    
    # We do use predictions of the doublet mode here for 16um and 32um
    # Can also use full mode predictions and rescale, but it will look a bit different
    ext <- ifelse(bs == 8, "", "_doubletmode")
    deconv_props <- readRDS(paste0("visium_hd_liver_combined/rds/liver",toupper(dataset), "_", bin_size_str, "_deconv_props", ext, ".rds"))
    
    plot_scatterbar(deconv_props, square_size,
                    visium_obj_roi =  visium_obj_roi,
                    color_palette = color_palette,
                    return_df = TRUE)
    
  }) %>% setNames(binsize_roi)
  
  # p_scatterbar <- purrr::map(p_props_list, "plot") 
  p_scatterbar <- lapply(binsize_roi, function(txt) {
    bs <- as.numeric(strsplit(txt, "_")[[1]][1])
    p_props_list[[txt]]$plot + ggtitle(paste0(bs, "\u00b5m"))
  }) %>% setNames(binsize_roi)
    
  props_df <- purrr::map(p_props_list, "df")
  
  # Make custom legend
  legend_square_size <- 0.9
  
  for (roi_i in 1:length(rois_oi)){
    spacing <- 0.25
    
    legend_df <- props_df[paste0(bin_sizes, "_", roi_i)] %>% 
      bind_rows(.id = "bin_size") %>% 
      ungroup() %>% 
      distinct(celltype, bin_size) %>% 

      mutate(bin_size = as.numeric(str_extract(bin_size, "([0-9]+)_", group=1))) %>% 
      mutate(x = case_when(bin_size == 8 ~ legend_square_size/6,
                           bin_size == 16 ~ legend_square_size/2,
                           bin_size == 32 ~ legend_square_size/2+legend_square_size/3),
             celltype = droplevels(factor(celltype, levels = celltype_order)),
             y = as.numeric(celltype) + (as.numeric(celltype)-1)*spacing, 
             celltype_label = celltype_labeller(celltype))
    
    legend_presence_df <- data.frame(
      x=c(legend_square_size/6, legend_square_size/2, legend_square_size/2+legend_square_size/3),
      y=c(-4, -3, -2),
      label=c("8\u00b5m", "16\u00b5m", "32\u00b5m")
    ) %>% mutate(y = y +(y+1)*spacing)
    
    p_legend_text_size <- 3
    p_legend <- ggplot(legend_df, aes(x = x, y = y)) +
      # Cell type legend
      geom_tile(aes(fill=celltype), height = legend_square_size, width = legend_square_size/3, show.legend = FALSE) +
      geom_text(data = legend_df %>% distinct(y, celltype_label), 
                aes(label=celltype_label, x=1.5), size=p_legend_text_size, hjust=0) +
      geom_text(data = data.frame(y=0, label="Cell types"), aes(label=label, y=y, x=0), size=p_legend_text_size, hjust=0) +
      # Presence legend
      geom_tile(data = legend_presence_df, aes(y=y), height=legend_square_size) +
      geom_tile(data = data.frame(y=legend_presence_df$y), aes(x=legend_square_size/2, y=y),
                height=legend_square_size, width=legend_square_size, fill=NA, color="black") +
      geom_text(data = legend_presence_df, aes(label=label, y=y, x=1.5), size=p_legend_text_size, hjust=0) +
      geom_text(data = data.frame(y=min(legend_presence_df$y)-1, label="Presence in ROI"),
                aes(label=label, y=y, x=0), size=p_legend_text_size, hjust=0) +
      scale_fill_manual(values = c(color_palette)) +
      xlim(0, 10) +
      theme_void() +
      scale_y_reverse() +
      coord_fixed()
    
    design <- "
    123
    145
    "
  
    row_1 <- p_roi + p_features_list[[paste0("Glul_", roi_i)]] + 
      p_features_list[[paste0("Hal_", roi_i)]] +
      p_features_list[[paste0("Spp1_", roi_i)]] +
      p_features_list[[paste0("Clec4f_", roi_i)]] +
      plot_layout(design=design, guides="collect") &
      theme(plot.title = element_text(hjust = 0.5, size=7),
            legend.key.size = unit(0.5, "cm"),
            legend.text = element_text(size=6),
            legend.title = element_text(size=6))
    row_1 <- wrap_elements(row_1)
    
    row_3 <- p_scatterbar[[paste0("16_", roi_i)]] + 
      p_scatterbar[[paste0("32_", roi_i)]] &
      theme(legend.position="none",
            plot.title = element_text(hjust = 0.5, size=7))
    
    row_23 <- wrap_elements(
      p_scatterbar[[paste0("8_", roi_i)]] / row_3 +
        plot_layout(heights=c(1.05, 0.5)) &
        theme(legend.position="none",
              plot.title = element_text(hjust = 0.5, size=7)))
    
    row_23_with_legend <- wrap_elements(row_23 + p_legend + 
                                          plot_layout(widths = c(0.75, 0.25)))
    
    
    all_rows <- row_1 / row_23_with_legend +
      plot_layout(heights = c(0.3, 0.7)) 
    
    ggsave(paste0(plot_path, "spatialbarplot_allres_featureplot", dataset, "_roi", roi_i, ".pdf"),
           all_rows, width = 8, height = 9, bg = "white")
    
  }
  
}