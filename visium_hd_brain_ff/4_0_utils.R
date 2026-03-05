#### Function to rotate image ####
library(cowplot)
library(grid)

rotate_image <- function(p, rot_angle) {
  gt <- ggplot_gtable(ggplot_build(p))
  panel_idx <- which(gt$layout$name == "panel")
  rot_vp <- viewport(angle = rot_angle)
  gt[["grobs"]][[panel_idx]] <- editGrob(gt[["grobs"]][[panel_idx]], vp = rot_vp)
  p_rot <- ggdraw() + draw_grob(gt)
  
  return(p_rot)
}

orig_xlim <- c(59.16006, 505.66585)
orig_ylim <- c(-60.19986, 386.12514)

# These numbers were determined from here
# orig_xlim <- ggplot_build(p_ncount)$layout$panel_scales_x[[1]]$range$range
# orig_ylim <- ggplot_build(p_ncount)$layout$panel_scales_y[[1]]$range$range
