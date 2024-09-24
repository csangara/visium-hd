# This script transfers region annotations to Visium slides

library(imager)
library(dplyr)
library(ggplot2)
library(httr)
library(jsonlite)

#### TRANSFERRING REGION ANNOTATIONS TO VISIUM DATA ####
# Load transformed atlas map (from QuickNII or VisuAlign)

coronal_visualign_map <- load.image("data/Visium_HD_MouseBrain/QUINT/tissue_8um_lowres_image_flat.png")

# Load Visium image
coronal_visium_img <- load.image("data/Visium_HD_MouseBrain/square_008um/spatial/tissue_lowres_image.png")

# Resize the atlas map as the Visium image size
dim(coronal_visium_img)[1:2]
dim(coronal_visualign_map)
coronal_visualign_map_resize <- resize(coronal_visualign_map, size_x = dim(coronal_visium_img)[1],
                                       size_y = dim(coronal_visium_img)[2])
plot(coronal_visualign_map)          # They seem similar
plot(coronal_visualign_map_resize)

# Load scale factors
scale_factors <- fromJSON("data/Visium_HD_MouseBrain/square_008um/spatial/scalefactors_json.json")

# Load Visium coordinates
coronal_visium_coords <- arrow::read_parquet("data/Visium_HD_MouseBrain/square_008um/spatial/tissue_positions.parquet") %>% 
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
rainbow_map <- jsonlite::fromJSON("data/Visium_HD_MouseBrain/QUINT/Rainbow 2017.json")
rainbow_map <- rainbow_map %>% mutate(rgb.val=rgb(red, green, blue, maxColorValue = 255))

# Number of regions in the dictionary vs in the picture
length(unique(rainbow_map$rgb.val))
length(unique(df$rgb.val))
unique(df$rgb.val) %in% unique(rainbow_map$rgb.val)

# Merge coordinates with the atlas map
coords_to_color <- merge(x = coronal_visium_coords %>%
                           mutate(x = round(pxl_col_in_lowres),
                                  y = round(pxl_row_in_lowres)),
                         y = df %>% select(x, y, rgb.val),
                         by = c("x", "y"))

# Merge that again with the region annotation
coords_to_region <- merge(x = coords_to_color,
                          y = rainbow_map %>% select(index, name, rgb.val),
                          by=c("rgb.val"))

# Quick and dirty: remove duplicated
coords_to_region <- coords_to_region[!duplicated(coords_to_region$barcode),]

# Let's check the regions interactively
library(plotly)
p <- ggplot(coords_to_region,aes(y=pxl_row_in_lowres, x=pxl_col_in_lowres, fill=name)) +
  geom_bin2d(bins=200) + scale_y_reverse() + coord_fixed(ratio=1) +
  theme_classic() +
  theme(legend.position="none")
p
ggplotly(p)

# Next step: convert the index to a broader hierarchy
library(igraph)

# Get entire structure graph
res <- GET("http://api.brain-map.org/api/v2/structure_graph_download/1.json")
content <- fromJSON(rawToChar(res$content))$msg

# Create starting graph of the root node
hierarchy <- make_graph()
hierarchy <- add_vertices(hierarchy, 1,
                          attr = list("id" = content$id,
                                      "st_level" = content$st_level,
                                      "acronym" = content$acronym,
                                      "parent_structure_id" = content$parent_structure_id,
                                      "graph_order" = content$graph_order))
  
to_traverse <- list(content$children[[1]])
while (length(to_traverse) > 0){
  current_child <- to_traverse[[1]]
  to_traverse <- to_traverse[-1]
  
  # Add child nodes
  hierarchy <- add_vertices(hierarchy, nrow(current_child),
                            attr = list("id" = current_child$id,
                                        "st_level" = current_child$st_level,
                                        "acronym" = current_child$acronym,
                                        "parent_structure_id" = current_child$parent_structure_id,
                                        "graph_order" = current_child$graph_order))
  
  # Add edges from each of the current child
  for (c in 1:nrow(current_child)){
    parent_vertex_id <- V(hierarchy)[V(hierarchy)$id == current_child[c,]$parent_structure_id]
    child_vertex_id <- V(hierarchy)[V(hierarchy)$id == current_child[c,]$id]
    hierarchy <- add_edges(hierarchy, edges = c(child_vertex_id, parent_vertex_id))
    
    # Add children
    if (!rlang::is_empty(current_child$children[[c]])){
      to_traverse <- c(to_traverse, current_child[c,]$children)
    }
  }
}

plot(hierarchy, layout = layout_as_tree(hierarchy, flip.y = FALSE))

# Now, we want to check different hierarchy labels to our annotations
all_regions <- (unique(coords_to_region$index)-1) %>% .[. >= 0]
all_lvls <- lapply(all_regions, function(ind) {
  # Get shortest path between current index and the root
  path <- shortest_paths(hierarchy, match(ind, V(hierarchy)$graph_order), 1)$vpath[[1]]
  
  # Find the parent whose st_level is from 3 to 8
  lapply(3:8, function(i) {
    lvl <- V(hierarchy)[V(hierarchy) %in% path & V(hierarchy)$st_level == i]
    
    if (rlang::is_empty(lvl)) {
      lvl <- list(acronym = NA, graph_order = NA)
      lvl_id <- NA
    } else {
      lvl_id <- as_ids(lvl)
    }
    return(data.frame(id = lvl_id, acronym = lvl$acronym, graph_order = lvl$graph_order))
  }) %>% bind_rows() %>%  set_rownames(paste0("lvl", 3:8)) %>% 
    tibble::rownames_to_column("type") %>% tidyr::pivot_wider(names_from = type, values_from = -type) %>% 
    mutate(index = ind+1, .before = 1) 

}) %>% bind_rows()

# Check number of unique regions in each level
length(unique(all_lvls$acronym_lvl6))
length(unique(all_lvls$acronym_lvl7))
length(unique(all_lvls$acronym_lvl8))

# Finally, we can append this information to the spot metadata!
coords_with_parent_regions <- merge(coords_to_region, all_lvls %>% select(index, acronym_lvl5:acronym_lvl8),
                                    by= "index")

ggplot(coords_with_parent_regions,aes(y=pxl_row_in_lowres, x=pxl_col_in_lowres, fill=acronym_lvl6)) +
  geom_bin2d(bins=500) + scale_y_reverse() + coord_fixed(ratio=1) +
  theme_classic() +
  theme(legend.position="none")
p
# ggplotly(p)

# Save the metadata
saveRDS(coords_with_parent_regions, "data/Visium_HD_MouseBrain/tissue_positions_with_annotations_008um.rds")
