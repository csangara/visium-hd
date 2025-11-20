library(Seurat)
library(tidyverse)

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


liver_atlas <- readRDS("~/spotless-benchmark/data/rds/liver_mouseStSt_noEC.rds")

atlas_df <- liver_atlas@meta.data %>% 
  # Remove spaces and dots from cell type names
  mutate(annot_vizgen = gsub(" ", "", annot)) %>%
  # Group cell types into those matching vizgen
  mutate(annot_vizgen = case_when(annot_vizgen %in% c("Fibroblasts") ~ "Stromalcells",
                                  grepl("Tcells|NKcells|Basophils|DC|Monocytes|ILC1s|Neutrophils|HsPCs", annot_vizgen) ~ "Otherimmunecells",
                                  # Others stay the same
                                  TRUE ~ annot_vizgen
  ))

digests <- c("exVivo", "inVivo", "nuclei")
proper_digest_names <- c("scRNA-seq\n(ex vivo digestion)", "scRNA-seq\n(in vivo digestion)", "snRNA-seq") %>%
  setNames(digests)

atlas_df %>% group_by(digest, annot_vizgen) %>%
  summarise(mean_counts = mean(nCount_RNA)) %>% 
  mutate(ratio_to_hepato = mean_counts[annot_vizgen == "Hepatocytes"]/mean_counts) %>% 
  #filter(!annot_vizgen %in% c("Hepatocytes", "Unknown")) %>% 
  filter(digest == "nuclei")
#group_by(digest) %>%

atlas_df %>%
  # Remove &nbsp; and dots from cell type names
  mutate(annot_cd45 = gsub("[^[:alnum:]]", "", annot_cd45)) %>%
  mutate(annot_cd45 = factor(annot_cd45, levels = rev(celltype_order))) %>% 
  ggplot(aes(y=annot_cd45, x=nCount_RNA, fill=annot_cd45)) +
  geom_boxplot(
    linewidth=0.15,
    outlier.size = 0.4,
    outlier.shape = 16,
    outlier.stroke = 0,
    show.legend = FALSE
  ) +
  stat_summary(fun=mean, geom="point", shape=23, size=0.5, stroke=0.25,
               color="black", fill="white") +
  facet_wrap(~digest,
             scales="free_x",
             labeller = as_labeller(proper_digest_names)
             ) +
  theme_minimal(base_size=7) +
  labs(x="Total UMIs per cell") +
  scale_fill_manual(values=color_palette,
                    name="Cell type") +
  scale_y_discrete(labels = celltype_labeller) +
  scale_x_continuous(labels = scales::label_number(scale_cut=scales::cut_short_scale())) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.25),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth=0.25),
        strip.text = element_text(size=7, face="bold"),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=6,margin=margin(t=5)))
ggsave("visium_hd_liver_combined/plots/boxplot_liver_atlas_counts_per_celltype.pdf",
       width=8, height=5, dpi=300)
