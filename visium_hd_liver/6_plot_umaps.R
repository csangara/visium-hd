library(Seurat)
library(tidyverse)
library(patchwork)

#### CREATE UMAP
nuclei_binned <- TRUE

color_palette_vishd <- c("Hepatocytes" = "#B4B5B5FF",
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

for (nuclei_binned in c(TRUE, FALSE)){
if (nuclei_binned){

  visium_hd_liver <- readRDS("data/Visium_HD_Liver/table_002um_stardist_segmentation_all_celltype_scores.rds")
  visium_hd_liver$cluster <- visium_hd_liver$clusters
  
  umap_meta <- bind_cols(visium_hd_liver@reductions$umap@cell.embeddings,
                         visium_hd_liver@meta.data[, c("annotation", "cluster"), drop=FALSE]) %>% 
    rename(UMAP1 = UMAP_1, UMAP2 = UMAP_2, celltype = annotation)
  
  # First run
  # Idents(visium_hd_liver) <- "clusters"
  # liver_markers <- FindAllMarkers(visium_hd_liver, only.pos = TRUE, min.pct = 0.1, )
  # saveRDS(liver_markers, "data/Visium_HD_Liver/table_002um_stardist_segmentation_cluster_markers.rds")
  
  liver_markers <- readRDS("data/Visium_HD_Liver/table_002um_stardist_segmentation_cluster_markers.rds")
  
} else {
  
  visium_hd_liver <- readRDS("data/Visium_HD_Liver/Visium_HD_Liver_008um.rds")
  removed_spots <- Cells(visium_hd_liver)[which(visium_hd_liver@meta.data[,"nCount_Spatial.008um"] < 100)]
  visium_hd_liver <- visium_hd_liver[, !Cells(visium_hd_liver) %in% removed_spots]
  deconv_props <- read.table("visium_hd_liver/Visium_HD_Liver_008um/proportions_rctd_Visium_HD_Liver_008um_annot_cd45",
                             header = TRUE)
  
  umap_meta <- read.csv("data/Visium_HD_Liver/Visium_HD_Liver_008um_umap_clusters.csv") %>% 
    mutate(cluster = as.character(cluster)) %>% 
    filter(!X %in% removed_spots) %>% 
    column_to_rownames("X")
  umap_meta$celltype <- colnames(deconv_props)[apply(deconv_props, 1, which.max)]
  
  rownames(deconv_props) <- Cells(visium_hd_liver)
  
  all(rownames(deconv_props) == rownames(umap_meta))
  all(Cells(visium_hd_liver) == rownames(umap_meta))
  
  visium_hd_liver <- NormalizeData(visium_hd_liver)
  visium_hd_liver <- AddMetaData(visium_hd_liver, umap_meta)
  visium_hd_liver$celltype_deconv <- colnames(deconv_props)[apply(deconv_props, 1, which.max)]
  
  # First run
  # Idents(visium_hd_liver) <- "cluster"
  # liver_markers <- FindAllMarkers(visium_hd_liver, only.pos = TRUE, min.pct = 0.1, )
  # saveRDS(liver_markers, "data/Visium_HD_Liver/Visium_HD_Liver_008um_cluster_markers.rds")
  
  liver_markers <- readRDS("data/Visium_HD_Liver/Visium_HD_Liver_008um_cluster_markers.rds")
}


umap_median <- umap_meta %>% 
  group_by(cluster) %>% 
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2))

p_umap_cluster <- ggplot(umap_meta) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = cluster), 
             size=0.2, stroke=0, shape=16) +
  geom_label(data = umap_median,
             aes(x = UMAP1, y = UMAP2, label = cluster, fill = cluster),
             color = "black", size=2) +
  scale_color_manual(values = RColorBrewer::brewer.pal(7, "Set2")) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(7, "Set2")) +
  guides(color = guide_legend(override.aes = list(size=1.5))) +
  ggtitle(paste0("Leiden clusters")) +
  theme_classic(base_size=7) +
  theme(plot.title = element_text(hjust = 0.5, size=7, face="bold"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

p_umap_cluster

p_celltype_title <- ifelse(nuclei_binned,
                     "Predicted cell type from marker scoring",
                     "Most abundant cell type from deconvolution")
p_celltype <- ggplot(umap_meta %>% group_by(celltype) %>% 
                       mutate(total_cells = n()) %>%
                       ungroup() %>% arrange(desc(total_cells))) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = celltype), 
             size=0.2, stroke=0, shape=16) +
  scale_color_manual(values = color_palette_vishd,
                     breaks=umap_meta$celltype %>% table %>%
                       sort(decreasing = TRUE) %>% .[1:10] %>% names) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  ggtitle(p_celltype_title) +
  theme_classic(base_size=7) +
  theme(plot.title = element_text(hjust = 0.5, size=7, face="bold"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

p_celltype

features_of_interest <- c("Glul", "Hal")
features_df <- data.frame(t(GetAssayData(visium_hd_liver)[features_of_interest, ])) %>% 
  rownames_to_column("cell_barcode") %>%
  inner_join(umap_meta %>% rownames_to_column("cell_barcode"),
             by="cell_barcode")

feature_plots <- lapply(features_of_interest, function(feat){
  umap_meta %>% rownames_to_column("spot") %>% 
    inner_join(GetAssayData(visium_hd_liver)[c(feat),,drop=FALSE] %>% t() %>% as.data.frame() %>%
                 rownames_to_column("spot"),
               by = "spot") %>%
    pivot_longer(cols = all_of(c(feat)), names_to = "gene", values_to = "expression") %>% 
    arrange(expression) %>%
    ggplot(aes(x=UMAP1, y=UMAP2, color=expression)) +
    geom_point(size=0.2, stroke=0, shape=16) +
    scale_color_gradientn(colors = c("lightgrey", "blue"), name="Expression") +
    ggtitle(feat) +
    theme_classic(base_size=7) +
    theme(plot.title = element_text(hjust = 0.5, size=7, face="bold"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          legend.direction = "horizontal",
          legend.position = "inside",
          legend.position.inside = c(0.1,0.95),
          legend.title = element_text(size=5, hjust=0.5),
          legend.title.position = "top",
          legend.text = element_text(size=5),
          legend.key.height = unit(0.3, "cm"),
          legend.key.width = unit(0.3, "cm"))
})

feature_plots[[1]]
feature_plots[[2]]

# DotPlot of top markers per cluster
nmarkers <- 5
nclusters <- visium_hd_liver$cluster %>% unique() %>% length()
top10_liver_markers <- liver_markers %>%
  group_by(cluster) %>%
  slice_head(n = nmarkers) %>%
  ungroup() %>% 
  mutate(cluster = factor(cluster, levels = 0:(nclusters-1))) %>% 
  arrange(cluster, desc(avg_log2FC))

dp <- DotPlot(visium_hd_liver, unique(top10_liver_markers$gene),
              group.by = "cluster") + coord_flip()

dotplot_df <- dp$data %>% 
  pivot_wider(names_from = id, names_prefix="g",
              id_cols = features.plot, values_from = c(avg.exp.scaled, pct.exp)) %>% 
  column_to_rownames("features.plot") %>%
  .[top10_liver_markers$gene,] %>% 
  mutate(gene_id = as.factor(row_number())) %>% 
  # pivot longer, but only separate groups
  pivot_longer(cols = -gene_id, names_to = c("measure", "cluster"),
               names_sep = "_", values_to = c("value")) %>% 
  pivot_wider(names_from = measure, values_from = value) %>% 
  mutate(gene_id_groups = rep(0:(nclusters-1), each=nclusters*nmarkers))

p_dotplot <- ggplot(dotplot_df) +
  geom_point(aes(x=cluster, y=forcats::fct_rev(gene_id), size=pct.exp, color=avg.exp.scaled)) +
  ggh4x::facet_wrap2(~gene_id_groups, ncol=1, scales="free_y", strip.position="left", dir="h") +
  scale_y_discrete(labels = function(x) top10_liver_markers$gene[as.numeric(x)]) +
  scale_size(range = c(0,2), name="Percent\nexpressed") +
  scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(9, "RdBu")), name = "Average\nexpression") +
  scale_x_discrete(name="Cluster", labels=0:(nclusters-1)) +
  theme_classic(base_size=7) +
  theme(strip.background = ggh4x::element_part_rect(side="r", linewidth=0.5, color="gray70"),
        strip.text.y.left = element_text(size=6, color="gray70", angle=0),
        strip.placement = "outside",
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=6),
        axis.text.y.left = element_text(size=5, margin=margin(l=5)),
        axis.title.x = element_text(size=6),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"))
p_dotplot

plots_list <- list(p_umap_cluster, p_celltype, feature_plots[[1]], feature_plots[[2]], p_dotplot)
# plots_8um <- plots_list
saveRDS(plots_list, ifelse(nuclei_binned,
                           "visium_hd_liver/rds/umap_plots_list_2um_stardist_binned.rds",
                           "visium_hd_liver/rds/umap_plots_list_8um.rds"))
}

plots_8um <- readRDS("visium_hd_liver/rds/umap_plots_list_8um.rds")
plots_2um <- readRDS("visium_hd_liver/rds/umap_plots_list_2um_stardist_binned.rds")

design <- "AABB
           AABB
           CCDD
           CCDD
           EEEE
           EEEE
           EEEE"

p_8um <- plots_8um[[1]] + plots_8um[[2]] + plots_8um[[3]] + plots_8um[[4]] + plots_8um[[5]] +
  plot_layout(design = design)

ggsave("visium_hd_liver/plots/umaps_8um.pdf",
       p_8um,
       width=6, height=8, dpi=300)

ggsave("visium_hd_liver/plots/umaps_8um.png",
       p_8um,
       width=6, height=8, dpi=300)


p_2um <- plots_2um[[1]] + plots_2um[[2]] + plots_2um[[3]] + plots_2um[[4]] + plots_2um[[5]] +
  plot_layout(design = design)


ggsave("visium_hd_liver/plots/umaps_2um_stardist_binned.pdf",
       p_2um,
       width=6, height=8, dpi=300)

ggsave("visium_hd_liver/plots/umaps_2um_stardist_binned.png",
       p_2um,
       width=6, height=8, dpi=300)

