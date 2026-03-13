library(synthesis)
library(sf)
library(spdep)
library(tidyverse)

# In this script, we want to see if Moran's I can be used to measure the 'spread' of cells
# This is because in segmented data (Vizgen), cells are clearly seen clustered together even
# when they have been artificially binned
# But in VisiumHD data, we see that marker genes are spread far apart

#### Helper functions ####
create_convex_blobs <- function(nobs=100, centers=5, sd=1, bbox=c(0,100), do.plot=TRUE, seed=123) {
  set.seed(seed)
  sim <- data.gen.blobs(
    nobs = nobs,
    features = 2,
    centers = centers,
    sd = sd,
    bbox = bbox,
    do.plot = do.plot
  )
  
  sim_points <- st_as_sf(as.data.frame(sim$x), coords = c("V1", "V2"))
  sim_points$cluster <- as.factor(sim$classes)
  
  sim_mbc <- sim_points %>%
    group_by(cluster) %>%
    summarize(geometry = st_combine(geometry)) %>%
    st_convex_hull()
  
  if (do.plot) {
    plot(sim_mbc, add=TRUE, border="red", lwd=2)
  }
  
  return(list(points=sim_points, convex_hulls=sim_mbc))
}

create_circles <- function(nobs=100, radius=5, bbox=c(0,100), do.plot=TRUE, seed=123) {
  # Pick nobs random points in a 2D space
  set.seed(seed)
  x <- runif(nobs, min=bbox[1], max=bbox[2])
  y <- runif(nobs, min=bbox[1], max=bbox[2])
  
  sim_points <- st_as_sf(data.frame(x=x, y=y), coords = c("x", "y")) %>% 
    # Turn to circle with st_buffer
    st_buffer(dist=radius) %>% 
    mutate(point_id=1:n())
  
  if (do.plot) {
    plot(sim_points)
  }
  return(sim_points)
}

create_grid_lw <- function(points, cellsize=c(8,8), bbox_range=c(0,100), do.plot=TRUE) {
  # Get bounding box of points
  bbox <- st_bbox(points)
  bbox["xmin"] <- min(bbox_range[1], bbox["xmin"])
  bbox["ymin"] <- min(bbox_range[1], bbox["ymin"])
  bbox["xmax"] <- max(bbox_range[2], bbox["xmax"])
  bbox["ymax"] <- max(bbox_range[2], bbox["ymax"])
  
  # Create a grid of specified cell size
  sim_grid <- st_make_grid(st_as_sfc(bbox), cellsize=cellsize) %>% 
    st_sf() %>% 
    mutate(grid_id=1:n())
  
  if (do.plot){
    plot(sim_grid)
    plot(points, add=TRUE)
  }
  
  # Create neighbors list from grid
  sim_nb <- poly2nb(sim_grid, queen=TRUE, snap=0.1)
  lw_sim <- nb2listw(sim_nb, style="W", zero.policy=TRUE)
  
  return(list(grid=sim_grid, lw=lw_sim, nb=sim_nb))
}

## Simulate Vizgen data - blobs = cells
# Create convex blobs  
vizgen_sim_points <- create_convex_blobs(nobs=1000, centers=10, sd=8, bbox=c(0,1000), do.plot=TRUE)$convex_hulls %>% 
  rename(point_id=cluster)
  
# If creating circle instead
# vizgen_sim_points <- create_circles(nobs=10, radius=8, bbox=c(0,1000), do.plot=TRUE)

# Create artificial grid
grid_lw <- create_grid_lw(vizgen_sim_points, cellsize=c(8,8), bbox_range=c(0,1000), do.plot=TRUE)

# Get overlapping areas of points and grid
vizgen_intersect <- st_join(grid_lw$grid,
                            vizgen_sim_points,
                            join=st_intersects)

vizgen_intersect_objects <- st_intersection(grid_lw$grid, vizgen_sim_points)

# Join with grid_id and point_id
vizgen_intersect_merge <- vizgen_intersect %>% st_set_geometry(NULL) %>% 
  left_join(vizgen_intersect_objects,
            by=c("grid_id", "point_id")) %>% 
  mutate(area = st_area(geometry),
         area_ratio = area / 64)

# Calculate Moran's I
vizgen_I <- moran.test(vizgen_intersect_merge$area_ratio,
                       grid_lw$lw, alternative="greater") 
vizgen_I #0.4

## Simulate VisiumHD data -> bigger points
vishd_sim_points <- create_convex_blobs(nobs=1000, centers=10, sd=32, bbox=c(0,1000), do.plot=TRUE)$convex_hulls %>% 
  rename(point_id=cluster)
vishd_grid_lw <- create_grid_lw(vishd_sim_points, cellsize=c(8,8), bbox_range=c(0,1000), do.plot=TRUE)

# Get overlapping areas of points and grid
vishd_intersect <- st_join(vishd_grid_lw$grid,
                            vishd_sim_points,
                            join=st_intersects)

vishd_intersect_objects <-  st_intersection(vishd_grid_lw$grid, vishd_sim_points)

# Join with grid_id and point_id
vishd_intersect_merge <- vishd_intersect %>% st_set_geometry(NULL) %>% 
  left_join(vishd_intersect_objects,
            by=c("grid_id", "point_id")) %>% 
  mutate(area = st_area(geometry)) %>% 
  group_by(grid_id) %>% 
  summarize(area = sum(area)) %>% 
  # Change area to 0.1 += gaussian noise if area >0, else area == 0
  mutate(area = ifelse(area > 0, 0.1 + rnorm(n(), mean=0.1, sd=0.05), 0))

vishd_I <- moran.test(vishd_intersect_merge$area,
                      vishd_grid_lw$lw, alternative="greater") 
vishd_I #0.8

# VisiumHD has a higher autocorrelation, because of the bigger area
# So Moran's I is not a good metric to measure gene 'spread'