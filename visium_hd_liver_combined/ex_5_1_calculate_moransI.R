source("visium_hd_liver_combined/0_utils.R")
library(spdep)
library(sf)

color_palette_combined <- c("Hepatocytes" = "#B4B5B5FF",
                            "Endothelialcells" = "#dca754",
                            "Cholangiocytes" = "#C61B84FF",
                            "Stromalcells" = "#79151e",
                            "Kupffercells" = "#5DA6DBFF",
                            "Otherimmunecells" = "#893A86FF",
                            "Bcells" = "#9C7EBAFF",
                            "Unknown" = "#191919",
                            "CentralVeinEndothelialcells" = "#FED8B1FF",
                            "PortalVeinEndothelialcells" = "#CC7722FF",
                            "Stellatecells" = "#A31A2AFF",
                            "Fibroblasts" = "#E45466FF")

#### MORAN'S I FOR CELL TYPES ####
binarize <- FALSE

# Vizgen Moran's I
vizgen <- read.csv(paste0("visium_hd_liver_combined/vizgen/vizgen_grid_",
                          bin_size, "um_celltype_proportions.csv")) %>% 
  # replace . in colname with nothing
  rename_with(~ gsub("\\.", "", .x)) %>% 
  # rename Kuppfercells to Kupffercells
  rename(Kupffercells = Kuppfercells)

grid_shp <- read_sf("visium_hd_liver_combined/vizgen/grid_8um.shp") %>% 
  mutate(FID = FID + 1) %>% 
  filter(FID %in% vizgen$grid_id)
grid_nb <- poly2nb(grid_shp, queen=TRUE, snap=1)
lw <- nb2listw(grid_nb, style="W", zero.policy=TRUE)

# Get Moran's I for all cell types
vizgen_moran_ct_df <- lapply(names(color_palette_vizgen), function(ct) {
  if (binarize) {
    vizgen_ct <- ifelse(vizgen[, ct] > 0, 1, 0)
  } else {
    vizgen_ct <- vizgen[, ct]
  }
  
  res <- moran.test(vizgen_ct, lw, alternative="greater") #p-value < 0.001
  
  stopifnot("all row names must be the same" = all(grid_shp$FID == vizgen$grid_id))
  
  data.frame(morans_I = res$estimate["Moran I statistic"], 
             expected_I = res$estimate["Expectation"], 
             variance_I = res$estimate["Variance"], 
             pvalue=res$p.value, dataset="Vizgen", celltype=ct
  )
}) %>% bind_rows() %>% `rownames<-`(NULL)


# VisiumHD Moran's I
vishd_shapes <- read_sf("visium_hd_liver_combined/rds/shapes/CAW009_square_008um.shp")
vishd_location <- read.csv("visium_hd_liver_combined/rds/shapes/CAW009_square_008um_locationid.csv", row.names=1) %>% 
  rownames_to_column("spot")

deconv_props_rank <- readRDS(paste0("visium_hd_liver_combined/rds/deconv_props_all.rds"))

deconv_props <- deconv_props_rank %>% 
  filter(dataset == "caw009", bin_size == "008um") %>% 
  ungroup() %>% 
  select(spot, celltype, proportion) %>% 
  tidyr::complete(spot, celltype, fill = list(proportion = 0))

vishd_shapes_sub <- vishd_shapes %>%
  inner_join(vishd_location %>% select(spot, location_id),
             by=c("location_i"="location_id")) %>% 
  filter(spot %in% unique(deconv_props$spot))
grid_vishd_nb <- poly2nb(vishd_shapes_sub, queen=TRUE, snap=0.1)
lw_vishd <- nb2listw(grid_vishd_nb, style="W", zero.policy=TRUE)

vishd_celltypes_oi <- c("Hepatocytes", "Cholangiocytes", "Bcells", "CentralVeinEndothelialcells",
                        "Kupffercells", "Fibroblasts", "PortalVeinEndothelialcells", "Stellatecells", "Tcells")

vishd_moran_ct_df <- lapply(vishd_celltypes_oi, function(ct) {
  
  deconv_props_ct <- deconv_props %>%
    filter(celltype == ct) %>% 
    # Arrange according to vishd_shapes_sub$spot
    arrange(match(spot, vishd_shapes_sub$spot))
  
  if (binarize) {
    deconv_props_ct$proportion <- ifelse(deconv_props_ct$proportion > 0, 1, 0)
  }
  
  stopifnot("all row names must be the same" = all(vishd_shapes_sub$spot == deconv_props_ct$spot))
  
  res <- moran.test(deconv_props_ct$proportion, lw_vishd, alternative="greater") #p-value < 0.001
  
  data.frame(morans_I = res$estimate["Moran I statistic"], 
             expected_I = res$estimate["Expectation"], 
             variance_I = res$estimate["Variance"], 
             pvalue=res$p.value, dataset="VisiumHD", celltype=ct)

}) %>% bind_rows() %>% `rownames<-`(NULL)

# Combine dataframes
moran_ct_df <- rbind(vishd_moran_ct_df, vizgen_moran_ct_df)

saveRDS(moran_ct_df, paste0("visium_hd_liver_combined/rds/moransI_celltypes",
                            ifelse(binarize, "_binarized", ""), ".rds"))

moran_ct_df_all <- lapply(c(TRUE, FALSE), function(b){
  readRDS(paste0("visium_hd_liver_combined/rds/moransI_celltypes",
                               ifelse(b, "_binarized", ""), ".rds")) %>% 
    mutate(binarized = b)
}) %>% bind_rows()

moran_ct_df_sub <- moran_ct_df_all %>% 
  # remove "Unknown", "Tcells", and "Otherimmunecells"
  filter(!celltype %in% c("Unknown", "Tcells", "Otherimmunecells")) %>%
  # Group 'fibroblasts', "Stellate cells", and "Stromal cells" together, keep the rest the same
  mutate(grouped_celltype = case_when(
    celltype %in% c("Fibroblasts", "Stellatecells", "Stromalcells") ~ "Fibro/Stromal/SC",
    celltype %in% c("Endothelialcells", "CentralVeinEndothelialcells", "PortalVeinEndothelialcells") ~ "Endothelialcells",
    TRUE ~ celltype
  )) %>% 
  mutate(grouped_celltype = factor(grouped_celltype, 
                                   levels=c("Hepatocytes", "Cholangiocytes", "Bcells", "Endothelialcells",
                                            "Kupffercells", "Fibro/Stromal/SC")),
         celltype = factor(celltype, levels=c("Hepatocytes", "Cholangiocytes", "Bcells", "CentralVeinEndothelialcells",
                                              "PortalVeinEndothelialcells", "Endothelialcells", "Kupffercells",
                                              "Fibroblasts", "Stellatecells", "Stromalcells"))) %>% 
  # Jitter y values only when the grouped_celltype contains multiple unique celltypes
  group_by(grouped_celltype) %>%
  mutate(y_axis = case_when(
    n_distinct(celltype) == 1 ~ as.numeric(grouped_celltype),
    TRUE ~ as.numeric(grouped_celltype) + seq(-0.2, 0.2, length.out=n_distinct(celltype))))

ggplot(moran_ct_df_sub,
       aes(y=y_axis, x=morans_I, fill=celltype, group=grouped_celltype, shape=dataset)) +
  geom_point(size=2.5, stroke=0.3) +
  theme_minimal(base_size=8) +
  labs(x="Moran's I", title=paste0("Spatial Autocorrelation of Cell types"),
       shape="Dataset") +
  scale_fill_manual(values=color_palette_combined, guide="none") +
  # scale shape: visiumhd = square, vizgen = triangle
  scale_shape_manual(values=c("VisiumHD"=22, "Vizgen"=21),
                     labels=c("Vizgen"="Gridded Vizgen")) +
  scale_y_continuous(breaks=sort(unique(moran_ct_df_sub$y_axis)),
                     labels=levels(moran_ct_df_sub$celltype),
                     # reverse
                     transform = "reverse") +
  scale_x_continuous(limits=c(0, 0.8)) +
  facet_wrap(~binarized, ncol=1,
             labeller = as_labeller(c(`TRUE` = "Binarized Proportions", `FALSE` = "Proportions"))
             ) +
  theme(legend.position="bottom",
        legend.title = element_text(size=6),
        legend.text = element_text(size=6, margin=margin(l=0)),
        legend.key.spacing.x = unit(0.3, "cm"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size=8, face="bold", hjust=0),
        panel.spacing.y = unit(0.5, "cm"),
        panel.border = ggh4x::element_part_rect(side="lb", fill=NA, linewidth=0.5),
        plot.title = element_text(hjust=0, size=8, face="bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=7),
        axis.ticks = element_line(color="black"))

ggsave(paste0("visium_hd_liver_combined/plots/moransI_celltypes_combined.pdf"),
       width=5, height=7)

#### CALCULATE LOCAL SPATIAL AUTOCORRELATION METRICS ####
# (Instead of the global metric)
# Try with KCs first
ct <- "Kupffercells"

# Vizgen data
vizgen <- read.csv(paste0("visium_hd_liver_combined/vizgen/vizgen_grid_",
                          bin_size, "um_celltype_proportions.csv")) %>% 
  # replace . in colname with nothing
  rename_with(~ gsub("\\.", "", .x)) %>% 
  # rename Kuppfercells to Kupffercells
  rename(Kupffercells = Kuppfercells)

grid_shp <- read_sf("visium_hd_liver_combined/vizgen/grid_8um.shp") %>% 
  mutate(FID = FID + 1) %>% 
  filter(FID %in% vizgen$grid_id)
grid_nb <- poly2nb(grid_shp, queen=TRUE, snap=1)
lw <- nb2listw(grid_nb, style="W", zero.policy=TRUE)
vizgen_ct <- vizgen[, ct]

# Local Moran's I
vizgen_localI <- localmoran(vizgen_ct, lw, alternative="greater")
hist(p.adjust(reslocal[,5], method="bonferroni"))

# Test out permutations + hotspot
vizgen_localIperm <- localmoran_perm(vizgen_ct, lw)
vizgen_Ih <- hotspot(vizgen_localIperm, Prname="Pr(z != E(Ii)) Sim", cutoff=0.05, p.adjust="none")
table(addNA(vizgen_Ih))

# Local Geary statistic
vizgen_localCperm <- localC_perm(vizgen_ct, lw)
vizgen_Ch <- hotspot(vizgen_localCperm, Prname="Pr(z != E(Ci)) Sim", cutoff=0.05, p.adjust="none")
table(addNA(vizgen_Ch))

# G local spatial statistics
vizgen_localGperm <- localG_perm(vizgen_ct, lw)
vizgen_Gh <- hotspot(vizgen_localGperm, Prname="Pr(z != E(Gi)) Sim", cutoff=0.05, p.adjust="none")
table(addNA(vizgen_Gh)) %>% prop.table

# VisiumHD data
vishd_shapes <- read_sf("visium_hd_liver_combined/rds/shapes/CAW009_square_008um.shp")
vishd_location <- read.csv("visium_hd_liver_combined/rds/shapes/CAW009_square_008um_locationid.csv", row.names=1) %>% 
  rownames_to_column("spot")
deconv_props_rank <- readRDS(paste0("visium_hd_liver_combined/rds/deconv_props_all.rds"))
deconv_props <- deconv_props_rank %>% 
  filter(dataset == "caw009", bin_size == "008um") %>% 
  ungroup() %>% 
  select(spot, celltype, proportion) %>% 
  tidyr::complete(spot, celltype, fill = list(proportion = 0))
vishd_shapes_sub <- vishd_shapes %>%
  inner_join(vishd_location %>% select(spot, location_id),
             by=c("location_i"="location_id")) %>% 
  filter(spot %in% unique(deconv_props$spot))
grid_vishd_nb <- poly2nb(vishd_shapes_sub, queen=TRUE, snap=0.1)
lw_vishd <- nb2listw(grid_vishd_nb, style="W", zero.policy=TRUE)

deconv_props_ct <- deconv_props %>%
  filter(celltype == ct) %>% arrange(match(spot, vishd_shapes_sub$spot)) %>% 
  pull(proportion)

# Same local statistics
vishd_localI <- localmoran(deconv_props_ct, lw_vishd, alternative="greater")
hist(p.adjust(vishd_localI[,5], method="bonferroni"))

vishd_localIperm <- localmoran_perm(deconv_props_ct, lw_vishd)
vishd_Ih <- hotspot(vishd_localIperm, Prname="Pr(z != E(Ii)) Sim", cutoff=0.05, p.adjust="none")
table(addNA(vishd_Ih))

vishd_localCperm <- localC_perm(deconv_props_ct, lw_vishd)
vishd_Ch <- hotspot(vishd_localCperm, Prname="Pr(z != E(Ci)) Sim", cutoff=0.05, p.adjust="none")
table(addNA(vishd_Ch))

vishd_localGperm <- localG_perm(deconv_props_ct, lw_vishd)
vishd_Gh <- hotspot(vishd_localGperm, Prname="Pr(z != E(Gi)) Sim", cutoff=0.05, p.adjust="none")
table(addNA(vishd_Gh))

# They don't work very well :(
ct_stats <- list(vizgen_I=vizgen_Ih, vizgen_C=vizgen_Ch, vizgen_G=vizgen_Gh,
     vishd_I=vishd_Ih, vishd_C=vishd_Ch, vishd_G=vishd_Gh)

# saveRDS(ct_stats, file="visium_hd_liver_combined/vizgen/local_stats_", ct, ".rds")
ct_stats <- readRDS(file="visium_hd_liver_combined/vizgen/local_stats_Kupffercells.rds")

for (i in 1:length(ct_stats)){
  print(names(ct_stats)[i])
  
  # Get letter for the type of stat (fifth letter fromt he end)
  stat_letter <- substr(names(ct_stats)[i], nchar(names(ct_stats)[i])-4, nchar(names(ct_stats)[i])-4)
  
  h <- hotspot(ct_stats[[i]], Prname=paste0("Pr(z != E(", stat_letter, "i)) Sim"), cutoff=0.05, p.adjust="none")
  print(table(addNA(h)))
}

# "vizgen_localIperm"
# Low-Low  Low-High High-High      <NA> 
#   1      6826      9890    221870 
# 
# "vizgen_localCperm"
# High-High        Low-Low Other Positive       Negative           <NA> 
#   8721              2           5307             21         224536 
#
# "vizgen_localGperm"
# Low   High   <NA> 
#   6822   9868 221897 
#
# "vishd_localIperm"
# Low-Low  High-Low  Low-High High-High      <NA> 
#   174        14     10216     30315    394012 
# 
# "vishd_localCperm"
# High-High        Low-Low Other Positive       Negative           <NA> 
#   32322            298           6918           2109         393084 
#
# "vishd_localGperm"
# Low   High   <NA> 
#   10526  30496 393709 


#### MORAN'S I for select genes ####
# Vizgen grid and neighbor list
grid_shp <- read_sf("visium_hd_liver_combined/vizgen/grid_8um.shp") %>% 
  mutate(FID = FID + 1)
grid_nb <- poly2nb(grid_shp, queen=TRUE, snap=1)
lw <- nb2listw(grid_nb, style="W", zero.policy=TRUE)

# VisiumHD grid and neighbor list
visium_obj <- readRDS("data/Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_008um.rds")
vishd_shapes <- read_sf("visium_hd_liver_combined/rds/shapes/CAW009_square_008um.shp")
vishd_location <- read.csv("visium_hd_liver_combined/rds/shapes/CAW009_square_008um_locationid.csv", row.names=1) %>% 
  rownames_to_column("spot")
vishd_shapes_sub <- vishd_shapes %>%
  inner_join(vishd_location %>% select(spot, location_id),
             by=c("location_i"="location_id"))
grid_vishd_nb <- poly2nb(vishd_shapes_sub, queen=TRUE, snap=0.1)
lw_vishd <- nb2listw(grid_vishd_nb, style="W", zero.policy=TRUE)

# Calculate Moran's I for KC marker genes
genes <- c("Sdc3", "Csf1r", "Kcnj16", "Clec4f")
gene <- "Sdc3"
morans_df <- lapply(genes, function(gene){
  
  if (gene != "Clec4f"){
    vizgen_gene <- read.csv(paste0("visium_hd_liver_combined/vizgen/vizgen_grid_",
                                   bin_size, "um_", gene, "_counts.csv"))
    
    vizgen_gene_filled <- vizgen_gene %>% right_join(grid_shp %>% st_drop_geometry() %>%
                                                       select(FID), by=c("index_right"="FID")) %>% 
      # Fill na with 0
      mutate(X0 = ifelse(is.na(X0), 0, X0))
    
    res_vizgen <- moran.test(vizgen_gene_filled$X0,lw, alternative="greater")
  } else {
    # Clec4f not in vizgen data
    res_vizgen <- list(estimate=c("Moran I statistic"=NA, "Expectation"=NA, "Variance"=NA), p.value=NA)
  }
  
  vishd_gene <- as.data.frame(GetAssayData(visium_obj, layer="counts")[gene, ]) %>%
    `colnames<-`(c("count")) %>% rownames_to_column("spot")
  
  # Check spot names match
  stopifnot("all row names must be the same" = all(vishd_shapes_sub$spot == vishd_gene$spot))
  
  res_vishd <- moran.test(vishd_gene$count,lw_vishd, alternative="greater") 
  
  data.frame(rbind(c(res_vishd$estimate, pvalue=res_vishd$p.value),
        c(res_vizgen$estimate, pvalue=res_vizgen$p.value)),
        dataset=c("VisiumHD", "Vizgen"),
        gene=gene)
}) %>% bind_rows()

# saveRDS(morans_df, file="visium_hd_liver_combined/vizgen/moransI_KC_genes.rds")
morans_df <- readRDS(file="visium_hd_liver_combined/vizgen/moransI_KC_genes.rds")

# Order gene by ascending difference between datasets
gene_order <- morans_df %>% group_by(gene) %>% 
  summarise(diff = abs(diff(Moran.I.statistic))) %>% 
  arrange(diff) %>% pull(gene)

p_moran_genes <- ggplot(morans_df %>% mutate(gene=factor(gene, levels=gene_order)),
       aes(y=gene, x=Moran.I.statistic, fill=dataset)) +
  geom_point(size=2, shape=21, stroke=0.3) +
  theme_minimal(base_size=7) +
  labs(x="Moran's I", title="Moran's I for KC Marker Genes", fill="Dataset") +
  scale_fill_manual(values=c("Vizgen"="#E69F00", "VisiumHD"="#56B4E9"),
                     labels=c("Vizgen"="Gridded Vizgen"))+
  theme(legend.title = element_text(size=6),
        legend.text = element_text(size=6, margin=margin(l=0)),
        legend.key.spacing.x = unit(0.3, "cm"),
        legend.position = "bottom",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"))
ggsave("visium_hd_liver_combined/plots/moransI_KC_genes.pdf",
       p_moran_genes,
       width=3, height=3) 

# Plot histogram of gene expressions
# To show that most of VisiumHD genes are zeroes
visium_obj <- readRDS("data/Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_008um.rds")
genes <- c("Sdc3", "Csf1r", "Kcnj16", "Clec4f")

vishd_genes_counts <- GetAssayData(visium_obj, layer="counts")[genes, ] %>% as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to="spot", values_to="count") 

vishd_genes_plots <- lapply(genes, function(g) {
  pg <- ggplot(vishd_genes_counts %>% filter(gene == g), aes(x=count)) +
    geom_histogram(bins=10, fill="#56B4E9") +
    theme_minimal(base_size=7) +
    labs(x="Counts", title=g) +
    ggbreak::scale_y_break(c(1e05, 3e05), scales="free") +
    # comma in y axis
    scale_y_continuous(labels=scales::comma) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=5)) +
    theme(panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          axis.ticks = element_line(color="black"),
          axis.title.y = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.line.y.right = element_blank())
  
  # ggsave(paste0("visium_hd_liver_combined/plots/histogram_", g, "_vishd.pdf"),
  #        pg, width=4, height=4, onefile=FALSE)
})


grid_shp <- read_sf("visium_hd_liver_combined/vizgen/grid_8um.shp") %>% 
  mutate(FID = FID + 1)
# Do the same for vizgen
vizgen_genes_plots <- lapply(genes[1:3], function(g) {
  vizgen_gene <- read.csv(paste0("visium_hd_liver_combined/vizgen/vizgen_grid_",
                                 bin_size, "um_", g, "_counts.csv"))
  
  vizgen_gene_filled <- vizgen_gene %>% right_join(grid_shp %>% st_drop_geometry() %>%
                                        select(FID), by=c("index_right"="FID")) %>% 
    # Fill na with 0
    mutate(X0 = ifelse(is.na(X0), 0, X0))
  
  pg <- ggplot(vizgen_gene_filled, aes(x=X0)) +
    geom_histogram(binwidth=1, fill="#E69F00") +
    theme_minimal(base_size=7) +
    labs(x="Counts", title=g) +
    ggbreak::scale_y_break(c(5e04, 2.5e05), scales="free") +
    # comma in y axis
    scale_y_continuous(labels=scales::comma) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=5)) +
    theme(panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          axis.ticks = element_line(color="black"),
          axis.title.y = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.line.y.right = element_blank())
  # ggsave(paste0("visium_hd_liver_combined/plots/histogram_", g, "_vizgen.pdf"),
  #        pg, width=4, height=4, onefile=FALSE)
})

legend <- cowplot::get_legend(p_moran_genes +
                                guides(fill = guide_legend(override.aes = list(shape = 22, size=4))) +
                                theme(legend.title.position = "top",
                                      legend.direction = "vertical",
                                      legend.key.spacing.y = unit(0, "cm")))

first_row <- cowplot::plot_grid(print(vishd_genes_plots[[1]]), print(vishd_genes_plots[[2]]),
                                print(vishd_genes_plots[[3]]), print(vishd_genes_plots[[4]]),
                                nrow=1)
second_row <- cowplot::plot_grid(print(vizgen_genes_plots[[1]]), print(vizgen_genes_plots[[2]]),
                                 print(vizgen_genes_plots[[3]]), legend,
                                 nrow=1)
cowplot::plot_grid(first_row, second_row, nrow=2)
ggsave(paste0("visium_hd_liver_combined/plots/histogram_KC_genes_visiumhd_vizgen.pdf"),
       width=8, height=6)

#### Calculate convex hull of connecting bins ####
# Instead of Moran's I, what if we calculate 'spread' of cell types based on
# the area of connecting bins?
# (Results of this section are probably not very trustworthy)

vishd_shapes <- read_sf("visium_hd_liver_combined/rds/shapes/CAW009_square_008um.shp")
vishd_location <- read.csv("visium_hd_liver_combined/rds/shapes/CAW009_square_008um_locationid.csv", row.names=1) %>% 
  rownames_to_column("spot")
deconv_props_rank <- readRDS(paste0("visium_hd_liver_combined/rds/deconv_props_all.rds"))
deconv_props <- deconv_props_rank %>% 
  filter(dataset == "caw009", bin_size == "008um") %>% 
  ungroup() %>% 
  select(spot, celltype, proportion) %>% 
  tidyr::complete(spot, celltype, fill = list(proportion = 0))

# Get connecting bins of cell types
# https://gis.stackexchange.com/questions/447556/grouping-polygons-based-on-clusters
for (ct in unique(deconv_props$celltype)){
  deconv_props_ct <- deconv_props %>%
    filter(celltype == ct, proportion > 0)
  vishd_shapes_sub_ct <- vishd_shapes %>%
    inner_join(vishd_location %>% filter(spot %in% deconv_props_ct$spot) %>% 
                 select(spot, location_id),
               by=c("location_i"="location_id"))
  
  # Buffer to make the corners touch, union to dissolve adjacent borders
  vishd_shapes_merge <- st_cast(st_union(st_buffer(vishd_shapes_sub_ct, 0.1)), "POLYGON")
  dissolved <- st_sf(vishd_shapes_merge) # Create a sf object from the geometries
  
  saveRDS(dissolved, paste0("visium_hd_liver_combined/rds/shapes/CAW009_square_008um_", ct, "_dissolved.rds"))
}

# Read scale factors
json_file <- "data/Visium_HD_Liver_CAW009/binned_outputs/square_008um/spatial/scalefactors_json.json"
json_data <- jsonlite::fromJSON(json_file)
vishd_scalefactor <- json_data$microns_per_pixel
# 1 pixel = 0.44223831 micron

vishd_dissolved <- lapply(unique(deconv_props$celltype), function(ct) {
  convex <- readRDS(paste0("visium_hd_liver_combined/rds/shapes/CAW009_square_008um_", ct, "_dissolved.rds")) %>% 
    st_convex_hull()
  data.frame(area = st_area(convex)*(vishd_scalefactor^2),
             celltype = ct)
}) %>% bind_rows()

ggplot(vishd_dissolved %>% filter(celltype != "Hepatocytes"),
       aes(x=area, y=celltype, color=celltype)) +
  scale_color_manual(values = color_palette) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "Area (µm<sup>2</sup>)") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = ggtext::element_markdown())

# Do the same for vizgen
vizgen <- read.csv(paste0("visium_hd_liver_combined/vizgen/vizgen_grid_",
                          bin_size, "um_celltype_proportions.csv")) %>% 
  # replace . in colname with nothing
  rename_with(~ gsub("\\.", "", .x)) %>% 
  # rename Kuppfercells to Kupffercells
  rename(Kupffercells = Kuppfercells)
grid_shp <- read_sf("visium_hd_liver_combined/vizgen/grid_8um.shp") %>% 
  mutate(FID = FID + 1) %>% 
  filter(FID %in% vizgen$grid_id)

vizgen_long <- vizgen %>%
  pivot_longer(-grid_id, names_to="celltype", values_to="proportion")

for (ct in colnames(vizgen)[-1]){
  grid_shp_ct <- grid_shp %>% 
    inner_join(vizgen_long %>% filter(celltype == ct, proportion > 0) %>%
                 select(grid_id),
               by=c("FID"="grid_id"))
  
  grid_shp_merge <- st_cast(st_union(st_buffer(grid_shp_ct, 0.1)), "POLYGON")
  dissolved <- st_sf(grid_shp_merge)
  
  saveRDS(dissolved, paste0("visium_hd_liver_combined/rds/shapes/vizgen_square_008um_", ct, "_dissolved.rds"))
}

vizgen_scalefactor <- 0.108
vizgen_dissolved <- lapply(colnames(vizgen)[-1], function(ct) {
  convex <- readRDS(paste0("visium_hd_liver_combined/rds/shapes/vizgen_square_008um_", ct, "_dissolved.rds")) %>% 
    st_convex_hull()
  data.frame(area = st_area(convex)*(vizgen_scalefactor^2),
             celltype = ct)
}) %>% bind_rows()

ggplot(vizgen_dissolved %>%
         filter(!celltype %in% c("Hepatocytes", "Unknown"), area < 50000),
       aes(x=area, y=celltype, color=celltype)) +
  scale_color_manual(values = color_palette_vizgen) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "Area (µm<sup>2</sup>)") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = ggtext::element_markdown())

# Original vizgen shapes
# Read in shape file
vizgen_ori_annot <- read.csv("visium_hd_liver_combined/vizgen/sdata_polyT_obs.csv") %>% 
  select(cells, annotation_own_score_genes) %>% 
  mutate(cells = as.character(cells)) %>% 
  rename(celltype = annotation_own_score_genes) %>%
  mutate(celltype = gsub(" ", "", celltype)) %>%
  mutate(celltype = ifelse(celltype == "Kuppfercells", "Kupffercells", celltype))

segmented_cells <- read_sf("visium_hd_liver_combined/vizgen/segmentation_mask_boundaries_roi.shp") %>% 
  left_join(vizgen_ori_annot %>% select(cells, celltype),
            by=c("index"="cells"))

# Get area of cells
segmented_cells_area <- segmented_cells %>%
  mutate(area = st_area(st_convex_hull(geometry))*vizgen_scalefactor^2) %>%
  data.frame()

ggplot(segmented_cells_area %>% filter(!celltype %in% c("Hepatocytes", "Unknown")),
       aes(x=area, y=celltype, color=celltype)) +
  scale_color_manual(values = color_palette_vizgen) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "Area (µm<sup>2</sup>)") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = ggtext::element_markdown())