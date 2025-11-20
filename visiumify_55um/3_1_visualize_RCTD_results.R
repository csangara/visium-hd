library(Seurat)
library(tidyverse)

dataset <- "_caw009" #or "_caw009" or ""
data_path <- paste0("data/Visium_HD_Liver", toupper(dataset), "/")
proportions_path <- paste0("visiumify_55um/Visium_HD_Liver", toupper(dataset))
plot_path <- paste0("visiumify_55um/plots/")

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

celltype_order <- names(color_palette)

celltype_labeller <- function(celltype_str){
  celltype_str %>% stringr::str_replace_all("Endothelialcells", "ECs") %>%
    stringr::str_replace_all("MonocytesMonocytederivedcells", "Mono & Mono-derived cells") %>%
    stringr::str_replace_all("([A-Za-z]+)(cells)$", "\\1 \\2") %>%
    stringr::str_replace_all("([a-z]{2,})([A-Z])", "\\1 \\2")
}

bin_size <- 56
bin_size_str <- sprintf("%03dum", bin_size)

visium_obj <- readRDS(paste0(data_path, "Visium_HD_Liver", toupper(dataset), "_",
                             bin_size_str, ".rds"))
dim(visium_obj) # 19059 genes x 9515 spots

ext <- "_annot_cd45"

deconv_props <- read.table(paste0(proportions_path, "_", bin_size_str,
                                  "/proportions_rctd_Visium_HD_Liver", toupper(dataset), "_",
                                  bin_size_str, ext),
                           header = TRUE)
dim(deconv_props) #9457 spots

removed_spots <- Cells(visium_obj)[which(visium_obj@meta.data[,paste0("nCount_Spatial.", bin_size_str)] < 100)]

# Check if removed rows + leftover rows == total rows (yes)
length(removed_spots) + dim(deconv_props)[1] == dim(visium_obj)[2]

# Subset visium_obj to only include spots that were not removed
visium_obj_subset <- visium_obj[, !(colnames(visium_obj) %in% removed_spots)]

dim(visium_obj_subset) # 19059 genes x 9457 spots

# Add rownames to deconv_props
rownames(deconv_props) <- colnames(visium_obj_subset)

deconv_props_df_all <- data.frame(deconv_props) %>% 
  pivot_longer(cols = everything(), names_to = "celltype",
               values_to = "proportion")

# Cell type proportions barplot
deconv_props_summ <- deconv_props_df_all %>% 
  group_by(celltype) %>% 
  summarise(agg_proportion = mean(proportion))

# Barplot ordered by total abundance
p <- ggplot(deconv_props_summ, aes(x=reorder(celltype, agg_proportion), y=agg_proportion)) +
  geom_bar(stat="identity") +
  # Add text of value, rounded to two digits, only if value is > 0
  geom_text(aes(label=ifelse(agg_proportion > 0.001, round(agg_proportion, 3), "")), nudge_y = 0.03, size=2) +
  coord_flip() +
  theme_minimal(base_size = 8) +
  theme(axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()) +
  ggtitle("Average cell type proportions across tissue")

p

# Load in data from spotless
datasets <- 1:4
props_visld <- lapply(datasets, function(ds) {
 
      read.table(paste0("~/spotless-benchmark/deconv_proportions/liver_mouseVisium_JB0", ds, "/proportions_rctd",
                        "_liver_mouseVisium_JB0", ds, "_noEC_annot_cd45"), header=TRUE, sep="\t") %>%
        # Still has . in colnames
        `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", ""))
}) %>% setNames(paste0("JB0", datasets)) %>%  reshape2::melt() %>% 
  `colnames<-`(c("celltype", "proportion", "slice"))

# Summarize mean proportions per slide
props_summ <- props_visld %>% group_by(celltype) %>%
  summarise(agg_proportion = mean(as.numeric(proportion))) %>% ungroup

# Barplot of deconv_props_summ and props_summ
bind_rows(deconv_props_summ %>% mutate(source = "VisiumHD_binned"),
          props_summ %>% mutate(source = "Visium")) %>%
  mutate(source = factor(source, levels = c("VisiumHD_binned", "Visium"))) %>% 
  ggplot(aes(y=reorder(celltype, agg_proportion), x=agg_proportion, fill=source)) +
  geom_bar(stat="identity", position = position_dodge()) +
  scale_y_discrete(labels = celltype_labeller) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 1, by=0.2)) +
  scale_fill_manual(values=c("Visium"="#ff7f0e", "VisiumHD_binned"="#1f77b4"),
                    limits=c( "Visium", "VisiumHD_binned"),
                    name="Dataset",
                    labels=c("VisiumHD_binned"="VisiumHD\n(binned to 56\u03bcm)",
                             "Visium"="Visium (55\u03bcm)")) +
  labs(x="Mean Proportion") +
  theme_bw(base_size = 8) +
  theme(panel.grid.major.y = element_blank(),
        panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.25),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth=0.25),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=6),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 5),
        legend.background = element_rect(fill = "white", color="black", linewidth=0.1),
        legend.position = "inside",
        legend.position.inside = c(0.85,0.1),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5))
ggsave(paste0(plot_path, "barplot_proportions_VisiumHD_binned_vs_Visium.pdf"),
       device = cairo_pdf,
       width=4.5, height=4, dpi=300)

# Nuc seq ground truth
# Filter out ABU21 (enriched for liver capsule)
# Also get coarser annotations
liver_nuc <- readRDS("~/spotless-benchmark/data/rds/liver_mouseStSt_snRNAseq_guilliams2022.rds")
liver_nuc <- liver_nuc@meta.data %>% 
  filter(sample != "ABU21") %>%  group_by(sample, annot_cd45) %>%
  count() %>% group_by(sample) %>% mutate(props=n/sum(n))

liver_nuc_summ <- liver_nuc %>% group_by(annot_cd45) %>% summarise(agg_proportion=mean(props)) %>% 
  # celltype, remove .&- replace with nothing
  mutate(celltype = stringr::str_replace_all(annot_cd45, "[ .&-]", ""))

vis_combined_df <- bind_rows(deconv_props_summ %>% mutate(source = "VisiumHD_binned"),
                             props_summ %>% mutate(source = "Visium"),
                             liver_nuc_summ %>% select(-annot_cd45) %>% mutate(source = "NucSeq"))

celltype_order <- vis_combined_df %>% group_by(celltype) %>% 
  summarise(summary = sum(agg_proportion)) %>%
  arrange(-summary) %>% 
  pull(celltype) %>% unique

vis_combined_agg <- vis_combined_df %>% 
  mutate(new_celltype = case_when(
    celltype %in% celltype_order[1:9] ~ celltype,
    TRUE ~ "Other")) %>% 
  group_by(new_celltype, source) %>%
  summarise(total_props = sum(agg_proportion)) %>% 
  mutate(new_celltype = factor(new_celltype, levels = c(celltype_order[1:9], "Other")))

celltype_format <- levels(vis_combined_agg$new_celltype) %>% 
  # Replace Endothelialcells with ECs
  stringr::str_replace_all("Endothelialcells", "ECs") %>%
  # Add space before the word "cells"
  stringr::str_replace_all("([A-Za-z]+)(cells)$", "\\1 \\2") %>% 
  # Add space before every capital letter that is not the first one
  stringr::str_replace_all("([a-z])([A-Z])", "\\1 \\2")

new_col <- c(color_palette %>% .[names(.) %in% celltype_order[1:9]],
             "Other" = "#191919")

ggplot(vis_combined_agg, aes(x=source, y=total_props, fill=new_celltype)) +
  geom_bar(stat="identity", width=0.75, position = position_dodge(width = 0.75)) +
  theme_minimal(base_size=8) +
  labs(x="Cell Type", y="Mean Proportion") +
  scale_fill_manual(values = new_col, name = "Celltype",
                    labels=celltype_format) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(labels=c("NucSeq" = "NucSeq (ground truth)", "Visium" = "Visium (55um)", "VisiumHD_binned" = "VisiumHD (binned to 56um)")) +
  guides(fill=guide_legend(nrow=2, title.position = "top", theme = theme(legend.byrow = TRUE))) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.25),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth=0.25),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(size=7, face="bold"),
        legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"),
        legend.key.spacing.x = unit(0.3, "cm"),
        legend.key.spacing.y = unit(0, "cm"),
        legend.title = element_text(size=6, hjust=0.5, margin=margin(b=5)),
        legend.text = element_text(size=6, margin = margin(l=3)))
ggsave(paste0(plot_path, "celltype_proportions_VisiumHD_binned_vs_Visium_vs_NucSeq.pdf"),
       width=7, height=3)

