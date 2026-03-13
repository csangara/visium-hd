library(ggalluvial)
library(tidyverse)

score_path <- "visium_hd_liver/score_genes"
gene_filters <- c("all", "q50", "q90")
score_type <- c("", "_iter")

scores_df <- lapply(gene_filters, function(gene_filter) {
  lapply(score_type, function(score) {
    path <- file.path(score_path, paste0("square_008um_atlas_", gene_filter, "_celltype_scores", score, "_props.csv"))
    read.csv(path) %>% rename(spot=X) %>% 
      as.data.frame() %>%
      mutate(gene_filter = gene_filter,
             score = ifelse(score == "", "raw", "iterative"))
  }) %>% bind_rows()
}) %>% bind_rows()

scores_wide <- scores_df %>% filter(score=="iterative") %>% select(-score) %>% 
  pivot_wider(names_from = "gene_filter", values_from = "annotation")

# Get all combinations of annotation shifts
scores_combi <-  dplyr::count(dplyr::group_by_all(dplyr::mutate_if(scores_wide %>% select(-spot), 
                is.numeric, function(x) factor(x, levels = as.character(sort(unique(x)))))), 
                name = "value")

ggplot(scores_combi,
       aes(y = value, axis1 = all, axis2 = q50, axis3 = q90)) +
  geom_alluvium() +
  geom_stratum(fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void()

is_alluvia_form(scores_combi, axes = 1:3, silent = TRUE)

scores_lodes <- to_lodes_form(scores_combi, axes = 1:3)
ggplot(scores_lodes,
       aes(y = value, alluvium=alluvium, x=x, stratum=stratum,
           fill=stratum)) +
  geom_flow() +
  geom_stratum() +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void()

library(wompwomp)
plot_alluvial(scores_wide, graphing_columns = c("all", "q50", "q90"), sorting_algorithm = "neighbornet",
              default_sorting = "decreasing", coloring_algorithm = "left",
              color_bands = FALSE, verbose = TRUE,
              min_text = 0.01, save_height = 24, save_width = 24, text_size = 35,
              axis_text_size = 25, axis_text_vjust = -5, text_width = 0.5)