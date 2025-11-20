# This script transfers region annotations to Visium slides

library(imager)
library(dplyr)
library(ggplot2)
library(httr)
library(jsonlite)

#### TRANSFERRING REGION ANNOTATIONS TO VISIUM DATA ####
# Load transformed atlas map (from QuickNII or VisuAlign)

coronal_visualign_map <- load.image("data/Visium_HD_MouseBrain_FF/QUINT/tissue_lowres_image_nl.png")

# Read parcellations file
parcellations_color <- read.csv("data/ABA_metadata/Allen-CCF-2020/20230630/views/parcellation_to_parcellation_term_membership_color.csv")
parcellations_acronym <- read.csv("data/ABA_metadata/Allen-CCF-2020/20230630/views/parcellation_to_parcellation_term_membership_acronym.csv")
parcellations <- inner_join(parcellations_color, parcellations_acronym,
                                 by="parcellation_index")
# Load Visium image
coronal_visium_img <- load.image("data/Visium_HD_MouseBrain_FF/binned_outputs/square_008um/spatial/tissue_lowres_image.png")

# Resize the atlas map as the Visium image size
dim(coronal_visium_img)[1:2]
dim(coronal_visualign_map)
coronal_visualign_map_resize <- resize(coronal_visualign_map, size_x = dim(coronal_visium_img)[1],
                                       size_y = dim(coronal_visium_img)[2])
plot(coronal_visualign_map)          # They seem similar
plot(coronal_visualign_map_resize)

# Load scale factors
scale_factors <- fromJSON("data/Visium_HD_MouseBrain_FF/binned_outputs/square_008um/spatial/scalefactors_json.json")

# Load Visium coordinates
coronal_visium_coords <- arrow::read_parquet("data/Visium_HD_MouseBrain_FF/binned_outputs/square_008um/spatial/tissue_positions.parquet") %>% 
  mutate(pxl_row_in_lowres = pxl_row_in_fullres * scale_factors$tissue_lowres_scalef,
         pxl_col_in_lowres = pxl_col_in_fullres * scale_factors$tissue_lowres_scalef) %>% 
  filter(in_tissue == 1)

ggplot(coronal_visium_coords, aes(y=pxl_row_in_lowres, x=pxl_col_in_lowres)) +
  geom_point(size=0.05) + coord_fixed(ratio=1) +
  scale_y_reverse()

# Plot coordinates on top of atlas map - looks ok
df <- as.data.frame(coronal_visualign_map_resize,wide="c") %>% mutate(rgb.val=rgb(c.1,c.2,c.3))
ggplot(df,aes(x,y)) + geom_raster(aes(fill=rgb.val)) +
  geom_point(inherit.aes=FALSE, data=coronal_visium_coords,
             aes(y=pxl_row_in_lowres, x=pxl_col_in_lowres), size=0.05, color="red", alpha=0.1) +
  scale_fill_identity() + scale_y_reverse() + coord_fixed(ratio=1)

# Load the "color to region name" dictionary file from https://www.nitrc.org/forum/message.php?msg_id=32830
#rainbow_map <- jsonlite::fromJSON("data/Visium_HD_MouseBrain/QUINT/Rainbow 2017.json")
# rainbow_map <- rainbow_map %>% mutate(rgb.val=rgb(red, green, blue, maxColorValue = 255))


# Number of regions in the dictionary vs in the picture
length(unique(parcellations_color$division_color))
length(unique(parcellations_color$structure_color))
length(unique(df$rgb.val))
unique(df$rgb.val) %in% unique(parcellations_color$structure_color)

# Which df$rgb.val are not in the parcellations?
setdiff(unique(df$rgb.val), unique(parcellations_color$structure_color))
"#66A83D" %in% unique(parcellations_color$substructure_color)

parcellations %>% filter(substructure_color == "#66A83D")

# Merge coordinates with the atlas map
coords_to_color <- merge(x = coronal_visium_coords %>%
                           mutate(x = round(pxl_col_in_lowres),
                                  y = round(pxl_row_in_lowres)),
                         y = df %>% select(x, y, rgb.val),
                         by = c("x", "y"))

# Check which colors are duplicated (fiber tracts and VS)
color_to_remove <- parcellations %>% distinct(structure_color, division) %>% 
  group_by(structure_color) %>%
  filter(n()>1) %>% pull(structure_color) %>% unique
parc_filtered <- parcellations %>% distinct(structure_color, division) %>%
  filter(!structure_color %in% color_to_remove) %>% 
  bind_rows(data.frame(structure_color = c(color_to_remove, "#66A83D"),
                       division = c(rep("fiber tracts/VS", 2), "HPF")))


# Merge that again with the region annotation
coords_to_region <- left_join(x = coords_to_color,
                          y = parc_filtered, 
                          by=c("rgb.val" = "structure_color"))

division_colors_df <- parcellations %>% distinct(division, division_color)
division_colors <- setNames(division_colors_df$division_color, division_colors_df$division)

# Let's check the regions interactively
library(plotly)
p <- ggplot(coords_to_region,aes(y=pxl_row_in_lowres, x=pxl_col_in_lowres, fill=division)) +
  geom_bin2d(bins=200) + scale_y_reverse() + coord_fixed(ratio=1) +
  scale_fill_manual(values = division_colors) +
  theme_classic() +
  theme()
p
ggplotly(p)

# Save the metadata
saveRDS(coords_to_region, "data/Visium_HD_MouseBrain_FF/tissue_positions_with_annotations_008um.rds")
