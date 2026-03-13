source("visium_hd_liver_combined/0_utils.R")

first_run <- TRUE

#### FIRST RUN - SAVE DECONV PROPS ####
# Combines proportions from both datasets at all three resolutions
# 8um is in doublet mode, 16 and 32um is in full mode
# But we also save 16 and 32um run in doublet mode
if (first_run){
for (dataset_name in c("sca002", "caw009")){
  dataset <- ifelse(dataset_name == "caw009", "_caw009", "")
  data_path <- paste0("data/Visium_HD_Liver", toupper(dataset), "/")
  proportions_path <- paste0("visium_hd_liver", dataset, "/Visium_HD_Liver", toupper(dataset))
  
  for (bin_size in bin_sizes){
    
    bin_size_str <- sprintf("%03dum", bin_size)
    
    visium_obj <- readRDS(paste0(data_path, "Visium_HD_Liver", toupper(dataset), "_",
                                 bin_size_str, ".rds"))
    
    for (ext in c("_annot_cd45", "_full")){
      if (ext == "_full" & bin_size == 8){
        next
      }
      
      # If _full mode, don't add extension
      file_ext <- ""
      if (bin_size %in% c(16, 32) & ext == "_annot_cd45"){ 
        file_ext <- "_doubletmode"
      }
      gc()
      deconv_props <- read.table(paste0(proportions_path, "_", bin_size_str,
                                        "/proportions_rctd_Visium_HD_Liver", toupper(dataset), "_",
                                        bin_size_str, ext),
                                 header = TRUE)
      
      removed_spots <- Cells(visium_obj)[which(visium_obj@meta.data[,paste0("nCount_Spatial.", bin_size_str)] < 100)]
      
      # Check if removed rows + leftover rows == total rows (yes)
      stopifnot(length(removed_spots) + dim(deconv_props)[1] == dim(visium_obj)[2])
      
      # Print number of removed spots, total spots, and percentage
      cat(paste(dataset_name, bin_size_str, "\n") )
      cat("Number of removed spots:", length(removed_spots), "\n")
      cat("Total spots:", dim(visium_obj)[2], "\n")
      cat("Percentage of removed spots:", 
          round(length(removed_spots) / dim(visium_obj)[2] * 100, 2), "%\n")
      
      # Subset visium_obj to only include spots that were not removed
      visium_obj_subset <- visium_obj[, !(colnames(visium_obj) %in% removed_spots)]
      
      # Add rownames to deconv_props
      rownames(deconv_props) <- colnames(visium_obj_subset)
      
      saveRDS(deconv_props, paste0("visium_hd_liver_combined/rds/", "liver", toupper(dataset), "_",
                                   bin_size_str, "_deconv_props", file_ext, ".rds"))
    }

  }
}
  
  # Combine all deconv props
  dataset_binsize <- paste0(rep(c("", "_CAW009"), each=length(bin_sizes)), "_", sprintf("%03dum", bin_sizes))
  deconv_props_all <- lapply(dataset_binsize, function(x) {
    deconv_props <- readRDS(paste0("visium_hd_liver_combined/rds/liver", x, "_deconv_props.rds"))

    data.frame(deconv_props) %>%
      rownames_to_column("spot") %>%
      pivot_longer(cols = -spot, names_to = "celltype",
                   values_to = "proportion") %>%
      mutate(dataset = ifelse(grepl("CAW009", x), "caw009", "sca002"),
             bin_size = gsub("_CAW009|_", "", x))
  }) %>% bind_rows

  deconv_props_rank <- deconv_props_all %>% filter(proportion > 0) %>%
    group_by(spot, dataset, bin_size) %>%
    mutate(rank = data.table::frankv(proportion, order=-1),
           doublet_type = case_when(n() == 1 ~ "singlet",
                                    n() > 1 ~ "doublet"))

  celltype_order <- deconv_props_rank %>%
    group_by(celltype, dataset, bin_size, rank, doublet_type) %>%
    summarise(total_counts = n()) %>%
    arrange(-total_counts) %>%
    pull(celltype) %>% unique

  deconv_props_rank <- deconv_props_rank %>%
    mutate(new_celltype = case_when(
      celltype %in% celltype_order[1:9] ~ celltype,
      TRUE ~ "Other")) %>%
    mutate(new_celltype = factor(new_celltype, levels = c("Other", rev(celltype_order[1:9]))))

  saveRDS(deconv_props_rank, paste0("visium_hd_liver_combined/rds/deconv_props_all.rds"))
  
  # Do this again for doublet mode only
  ds_bs_mode <- paste0(dataset_binsize,  "_deconv_props", rep(c("", "_doubletmode", "_doubletmode"), 2))
  
  deconv_props_all <- lapply(ds_bs_mode, function(x) {
    deconv_props <- readRDS(paste0("visium_hd_liver_combined/rds/liver", x, ".rds"))
    
    data.frame(deconv_props) %>%
      rownames_to_column("spot") %>%
      pivot_longer(cols = -spot, names_to = "celltype",
                   values_to = "proportion") %>%
      mutate(dataset = ifelse(grepl("CAW009", x), "caw009", "sca002"),
             bin_size = str_extract(x, "[0-9]+um"))
  }) %>% bind_rows
  
  deconv_props_rank <- deconv_props_all %>% filter(proportion > 0) %>%
    group_by(spot, dataset, bin_size) %>%
    mutate(rank = data.table::frankv(proportion, order=-1),
           doublet_type = case_when(n() == 1 ~ "singlet",
                                    n() > 1 ~ "doublet"))
  
  celltype_order <- deconv_props_rank %>%
    group_by(celltype, dataset, bin_size, rank, doublet_type) %>%
    summarise(total_counts = n()) %>%
    arrange(-total_counts) %>%
    pull(celltype) %>% unique
  
  deconv_props_rank <- deconv_props_rank %>%
    mutate(new_celltype = case_when(
      celltype %in% celltype_order[1:9] ~ celltype,
      TRUE ~ "Other")) %>%
    mutate(new_celltype = factor(new_celltype, levels = c("Other", rev(celltype_order[1:9]))))
  
  saveRDS(deconv_props_rank, paste0("visium_hd_liver_combined/rds/deconv_props_all_doubletmode.rds"))
}


#### PRINT TABLE OF MEAN PROPORTIONS ####
deconv_props_rank <- readRDS(paste0("visium_hd_liver_combined/rds/deconv_props_all.rds"))

dataset_binsize <- deconv_props_rank %>% ungroup() %>% 
  distinct(dataset, bin_size) %>% unite(dataset_bin, sep = "_") %>% 
  pull(dataset_bin)

mean_props <- lapply(dataset_binsize, function(x) {
  ds <- str_split(x, "_")[[1]][1]
  bin <- str_split(x, "_")[[1]][2]
  deconv_props_rank %>% ungroup() %>% 
    filter(dataset == ds, bin_size == bin) %>% 
    select(spot, celltype, proportion) %>% 
    pivot_wider(names_from = celltype, values_from = proportion, values_fill = 0) %>% 
    pivot_longer(cols = -spot, names_to = "celltype", values_to = "proportion") %>% 
    group_by(celltype) %>%
    summarise(mean_proportion = mean(proportion)) %>% 
    mutate(dataset = ds,
           bin_size = bin)
}) %>% bind_rows() %>% 
  tidyr::complete(dataset, bin_size, celltype,
                  fill = list(mean_proportion = 0))

mean_props %>% arrange(desc(dataset), bin_size) %>%
  pivot_wider(names_from = c(dataset, bin_size), values_from = mean_proportion) %>% 
  # get 2 significant digits in scientific notation
  mutate(across(where(is.numeric), ~ signif(., digits = 2))) %>%
  # Replace Endothelialcells with ECs in celltype
  mutate(celltype = stringr::str_replace_all(celltype, "Endothelialcells", "ECs")) %>% 
  # Add space before the word "cells"
  mutate(celltype = stringr::str_replace_all(celltype, "([A-Za-z]+)(cells)$", "\\1 \\2")) %>%
  write_csv(paste0("visium_hd_liver_combined/tables/deconv_props_mean.csv"))

mean_props_vizgen <- mean_props %>%
  # Remove spaces and dots from cell type names
  mutate(annot_vizgen = gsub(" ", "", celltype)) %>%
  group_annot_to_vizgen(annot_vizgen, "annot_vizgen") %>% 
  mutate(annot_vizgen = factor(annot_vizgen, levels = names(color_palette_vizgen)),
         dataset = factor(dataset, levels = c("sca002", "caw009"), labels = c("SCA002", "CAW009"))) %>% 
  # summarise to get mean proportion per annot_vizgen
  group_by(dataset, bin_size, annot_vizgen) %>%
  summarise(mean_proportion = sum(mean_proportion))

## PROPORTION BARPLOT - Figure 4.4d ##
p_mean_props <- ggplot(mean_props_vizgen %>% filter(bin_size == "008um"),
                       aes(x=dataset, y=mean_proportion, fill=annot_vizgen)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=color_palette_vizgen, name="Cell type") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 1, 0.2)) +
  guides(fill=guide_legend(nrow=1, title.position = "top")) +
  labs(title = "Mean cell type proportions in 8\u03bcm VisiumHD data",
       x="Cell Type", y="Mean Proportion") +
  theme_minimal(base_size=8) +
  theme_barplot +
  theme(axis.text.x = element_text(size=7),
        legend.position = "none")

design <- "####
           ####
           ##EE"

p_mean_props + plot_layout(design = design)
ggsave(paste0(plot_path, "vizgen_grid_with_props_UMIs_d.pdf"),
       device = cairo_pdf,
       width=8, height=6, dpi=300)

mean_props_factor <- mean_props %>% 
  mutate(dataset = factor(dataset, levels = c("sca002", "caw009"), labels = c("SCA002", "CAW009")),
         celltype = factor(celltype, levels = rev(names(color_palette))),
         bin_size = factor(bin_size, levels = c("008um", "016um", "032um"),
                           labels = c("8\u03bcm", "16\u03bcm", "32\u03bcm")))

## SUPPLEMENTARY FIGURE 4.4 - MEAN PROPS OF ALL CELL TYPES ##
ggplot(mean_props_factor,
       aes(y=celltype, x=mean_proportion, fill=celltype)) +
  geom_bar(stat="identity") +
  facet_grid(bin_size~dataset) +
  labs(y="Cell Type", x="Mean Proportion") +
  scale_fill_manual(values=color_palette, name="Cell type") +
  scale_y_discrete(labels = celltype_labeller) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size=8) +
  theme_barplot_facet +
  theme(strip.background.y = ggh4x::element_part_rect(side = "l", fill=NA, linewidth=0.5),
        strip.text.y.right = element_text(angle=0),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=5),
        legend.position = "none")

ggsave(paste0(plot_path, "barplot_deconv_mean_props_all_celltypes.pdf"),
       device = cairo_pdf,
       width=6, height=6, dpi=300)

## Calculate other statistics ##
# Calculate mean per cell type, dataset, bin size
deconv_props_rank %>% 
  group_by(celltype, dataset, bin_size) %>% 
  summarise(mean_proportion = mean(proportion)) %>% 
  arrange(celltype, dataset, bin_size) %>% 
  print(n=Inf)

# How many percent of spots are hepatocytes present?
deconv_props_rank %>% 
  filter(doublet_type == "doublet") %>% 
  group_by(dataset, bin_size) %>% 
  summarise(total_spots = n_distinct(spot),
            total_hepatocytes = sum(celltype == "Hepatocytes")) %>% 
  mutate(percent_hepatocytes = total_hepatocytes / total_spots * 100)


