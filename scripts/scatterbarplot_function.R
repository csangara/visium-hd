library(dplyr)
library(ggplot2)

#' Creates a scatterbarplot of a region of interest
#' 
#' The predicted proportions will be subset according to the observations present in `visium_obj_roi`.
#' 
#' @param deconv_props Dataframe of proportions, with observations in the rows and cell types in the column. Will be subset based on rownames to only contain observations present in `visium_obj_roi`.
#' @param square_size Square size of each bin. Can normally be obtained from `visium_obj_roi@images[["slice1.008um"]]@scale.factors$spot`.
#' @param ncelltypes The number of cell types to be plotted per bin. If `Inf`, plot all predicted cell types. As this may return a crowded plot, it is recommended that max. 2 cell types be plotted per bin. Proportions will be rescaled accordingly such that the total proportions of each spot sums to 1.
#' @param visium_obj_roi Seurat object that has been subset to region of interest. If NULL, will assume that `deconv_props` is in the correct format (see Details)
#' @param return_df Whether or not to also return the dataframe that is used for plotting.
#' @param color_palette Color palette to be used for plotting.
#' @returns A ggplot object, or also the underlying dataframe (if `return_df=TRUE`)
#' @details If `visium_obj_roi` is NULL, `deconv_props` must be in a long format, that is, each row consisting of the following columns: spot, celltype, proportion, x, y. The dataframe should also be filtered to the region of interest, and only include proportion > 0.

plot_scatterbar <- function(deconv_props,
                            square_size,
                            ncelltypes=Inf,
                            visium_obj_roi=NULL,
                            return_df=FALSE,
                            color_palette=NULL){
  
  stopifnot("ncelltypes must be 1 or higher." = ncelltypes >= 1)
  
  if (!is.null(visium_obj_roi)){
    # Filter deconv proportions to only present cells
    deconv_props_df <- deconv_props %>% 
      tibble::rownames_to_column("spot") %>%
      tidyr::pivot_longer(cols = -spot, names_to = "celltype",
                          values_to = "proportion") %>% 
      dplyr::filter(proportion > 0) %>% 
      # Add coordinates
      dplyr::right_join(Seurat::GetTissueCoordinates(visium_obj_roi),
                        by = c("spot" = "cell"))
  } else {
    deconv_props_df <- deconv_props
  }
  
  if (!is.infinite(ncelltypes)){
    cat(paste0("Limiting predictions to the top ", ncelltypes, " cell types and rescaling the proportions..."))
    
    deconv_props_df <- deconv_props_df %>% 
      dplyr::group_by(spot) %>%
      # Arrange each spot by descending proportion
      dplyr::arrange(spot, dplyr::desc(proportion)) %>% 
      dplyr::filter(dplyr::row_number() <= ncelltypes) %>% 
      # Rescale
      dplyr::mutate(proportion = proportion/sum(proportion))
  }

  # Create barplot
  deconv_props_df_barplot <- deconv_props_df %>% 
    dplyr::group_by(spot) %>%
    # Arrange each spot by descending proportion
    dplyr::arrange(spot, dplyr::desc(proportion)) %>% 
    dplyr::mutate(group = paste0(spot, "_", dplyr::row_number()),
                  cumu_prop = cumsum(proportion)) %>% 
    # The coords from GetTissueCoordinates are the center of the square
    # So we want to draw four corners of the squares instead (x1, y1) until (x4, y4)
    dplyr::mutate(x1 = y - (square_size/2) + (cumu_prop*square_size),
                  y1 = x - square_size / 2,
                  x2 = y - (square_size/2) + (cumu_prop*square_size),
                  y2 = x + square_size / 2,
                  x3 = y - (square_size/2) + ((cumu_prop - proportion)*square_size),
                  y3 = x + square_size / 2,
                  x4 = y - (square_size/2) + ((cumu_prop - proportion)*square_size),
                  y4 = x - square_size / 2) %>% 
    dplyr::rename(coord_x = x, coord_y = y) %>%
    tidyr::pivot_longer(cols = c(x1, y1, x2, y2, x3, y3, x4, y4),
                        names_to = c(".value", "corner"),
                        names_pattern = "(x|y)([1-4])")
  
  p <- ggplot(deconv_props_df_barplot, aes(x = x, y = y)) +
    # Main barplot
    geom_polygon(aes(fill = celltype, group = group), show.legend = TRUE) +
    # White borders
    geom_tile(data = deconv_props_df %>% distinct(x, y),
              aes(x = y, y = x), height = square_size, width = square_size,
              fill = NA, color = "white", inherit.aes = FALSE) +
    theme_void(base_size = 7) +
    scale_y_reverse() +
    coord_fixed() +
    guides(fill = guide_legend(ncol=1)) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=5),
          legend.key.size = unit(0.3, "cm"))
  
  if (!is.null(color_palette)){
    p <- p + scale_fill_manual(values = color_palette)
  }
  
  if (return_df){
    return(list(plot=p, df=deconv_props_df_barplot))
  } else{
    return(p)
  }
}
