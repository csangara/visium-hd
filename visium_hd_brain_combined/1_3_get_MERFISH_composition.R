source("visium_hd_brain_combined/0_brain_init.R")

first_run <- FALSE
if (first_run){
  merfish_zeng <- read.csv("data/ABA_metadata/MERFISH-C57BL6J-638850-CCF/20231215/views/cell_metadata_with_parcellation_annotation.csv")
  saveRDS(merfish_zeng, "visium_hd_brain_combined/rds/merfish_zeng.rds")
  
  parcellations_color <- read.csv("data/ABA_metadata/Allen-CCF-2020/20230630/views/parcellation_to_parcellation_term_membership_color.csv")
  parcellations_acronym <- read.csv("data/ABA_metadata/Allen-CCF-2020/20230630/views/parcellation_to_parcellation_term_membership_acronym.csv")
  parcellations <- inner_join(parcellations_color, parcellations_acronym,
                              by="parcellation_index")
  
  for (dsi in 1:2){
    merfish_zhuang <- read.csv(paste0("data/ABA_metadata/Zhuang-ABCA-", dsi, "/20241115/views/cell_metadata_with_cluster_annotation.csv"))
    merfish_zhuang_ccf <- read.csv(paste0("data/ABA_metadata/Zhuang-ABCA-", dsi, "-CCF/20230830/ccf_coordinates.csv"))
    merfish_zhuang_with_parcellation <- merfish_zhuang %>%
      inner_join(merfish_zhuang_ccf %>% rename(x_ccf = x, y_ccf = y, z_ccf = z), by="cell_label") %>%
      inner_join(parcellations, by="parcellation_index")
    
    saveRDS(merfish_zhuang_with_parcellation, paste0("visium_hd_brain_combined/rds/merfish_zhuang", dsi, "_with_parcellation.rds"))
            
  }
  
  rm(merfish_zhuang, merfish_zhuang_ccf, merfish_zhuang_with_parcellation); gc()
  
  merfish_zhuang1 <- readRDS("visium_hd_brain_combined/rds/merfish_zhuang1_with_parcellation.rds")
  merfish_zhuang2 <- readRDS("visium_hd_brain_combined/rds/merfish_zhuang2_with_parcellation.rds")
  
  # Combined dataframe
  merfish_zeng <- merfish_zeng %>%
    rename_with(~str_replace_all(., "^parcellation_", "")) %>% 
    rename(parcellation_index = index) %>% 
    rename_with(~str_replace_all(., "_section$", ""))
  
  common_cols <- intersect(colnames(merfish_zeng), colnames(merfish_zhuang1))
  
  merfish_cells <- bind_rows(
    merfish_zeng %>% select(all_of(common_cols)) %>%
      mutate(source = "Zeng"),
    merfish_zhuang1 %>% select(all_of(common_cols)) %>%
      mutate(source = "Zhuang1",
             cell_label = as.character(cell_label),
             donor_sex = "F"),
    merfish_zhuang2 %>% select(all_of(common_cols)) %>%
      mutate(source = "Zhuang2",
             cell_label = as.character(cell_label))
  )
  
  saveRDS(merfish_cells, "visium_hd_brain_combined/rds/merfish_cells_3datasets_combined.rds")
  
} else{
  merfish_cells <- readRDS("visium_hd_brain_combined/rds/merfish_cells_3datasets_combined.rds")
  
}

# How many sections per dataset?
merfish_cells %>% distinct(source, brain_section_label) %>% 
  count(source)

merfish_cells <- merfish_cells %>% 
  # Only get the last numbers after the dot
  mutate(brain_section_id = str_replace_all(brain_section_label, ".*\\.(\\d+)$", "\\1")) %>% 
  mutate(class = str_replace_all(class, "[- ]", ""))

library(ggridges)
ggplot(merfish_cells, aes(y=brain_section_id, x=x_ccf)) +
  geom_density_ridges() +
  facet_wrap(~source, scales = "free_y") +
  theme_bw()

# Plot distribution over all brain sections
# Check min, max, and difference of x_ccf
merfish_cells %>% 
  group_by(source, brain_section_id) %>% 
  mutate(mean_x_ccf = mean(x_ccf)) %>%  ungroup() %>% 
  summarise(min_x = min(mean_x_ccf),
            max_x = max(mean_x_ccf))

merfish_cells_bins <- merfish_cells %>% 
  group_by(source, brain_section_id) %>% 
  mutate(mean_x_ccf = mean(x_ccf)) %>%
  mutate(x_ccf_bins = cut(mean_x_ccf, breaks=seq(0.3, 12.6, by=0.2)))

merfish_cells_summ <- merfish_cells_bins %>%
  group_by(class, source, division, x_ccf_bins) %>% 
  summarise(n = n()) %>% 
  group_by(division, source, x_ccf_bins) %>% 
  mutate(prop = n / sum(n),
         class = factor(class, levels = names(celltype_colors)),
         x_ccf_bins = forcats::fct_rev(x_ccf_bins))

# How many brain sections are in a bin?
merfish_cells_bins %>% 
  group_by(source, brain_section_id) %>% 
  distinct(source, brain_section_id, x_ccf_bins) %>% 
  group_by(source, x_ccf_bins) %>% 
  summarise(n_sections = n()) %>% 
  group_by(source) %>%
  summarise(mean_n_sections = mean(n_sections),
            min_n_sections = min(n_sections),
            max_n_sections = max(n_sections))

merfish_cells %>% distinct(source, brain_section_label) %>% 
  count(source)

region_annotations <- readRDS(paste0("data/Visium_HD_MouseBrain_FF/tissue_positions_with_annotations_008um.rds"))
common_regions <- intersect(unique(region_annotations$division),
                            merfish_cells$division) %>% 
  # remove any that is 'unassigned'
  .[!grepl("unassigned", .)]

levels(merfish_cells_summ$x_ccf_bins)
ggplot(merfish_cells_summ %>%
         filter(division %in% common_regions),
         aes(x = prop, y = x_ccf_bins, fill = class)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE)) +
  facet_grid(source ~ division, drop=FALSE, labeller = labeller(division=proper_region_names)) +
  scale_fill_manual(values = celltype_colors, labels=proper_celltype_names) +
  scale_y_discrete(breaks=paste0("(", seq(0.9, 9.9, by=1), ",", seq(1.1, 10.1, by=1), "]"),
                   labels = 1:10,
                   limits = levels(merfish_cells_summ$x_ccf_bins)[11:61]) +
  scale_x_continuous(breaks=seq(0, 1, by=0.2), limits=c(0,1), labels=c(0, "", 0.4, "", 0.8, 1),
                     expand = expansion(mult = c(0, 0.02))) +
  theme_bw(base_size=8) +
  guides(fill=guide_legend(nrow=3)) +
  labs(y = "Allen Mouse Brain CCF x-axis", fill = "Cell type") +
  theme(legend.position="bottom",
        strip.background = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = ggh4x::element_part_rect(side = "l", fill = NA, linewidth = 0.5),
        panel.spacing.x = unit(5, "mm"),
        panel.spacing.y = unit(5, "mm"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=6),
        axis.text = element_text(size=5),
        legend.key.size = unit(0.3, "cm"),
        legend.key.spacing.x = unit(0.15, "cm"),
        legend.key.spacing.y = unit(0, "cm"),
        #legend.title = element_text(size=6, hjust=0.5, margin=margin(b=5)),
        legend.title = element_text(size=6, hjust=0.5, angle=90),
        legend.text = element_text(size=5, margin = margin(l=3)))

ggsave("visium_hd_brain_combined/plots/barplot_merfish_celltype_by_xccf_bins_by_region.pdf",
       device = cairo_pdf,
       width = 9, height = 6, dpi=300)

ggplot(merfish_cells_summ, aes(x = mean_x_ccf, y=source)) +
  geom_point() +
  theme_minimal()


