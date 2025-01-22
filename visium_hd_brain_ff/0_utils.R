library(RColorBrewer)
library(cowplot)
library(grid)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

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
