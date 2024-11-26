library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(caret)
library(precrec)

# First, let's read in scrublet and RCTD results, and compare them
doublet_props <- read.table(paste0("visium_hd_brain/Visium_HD_MouseBrain_008um/proportions_rctd_Visium_HD_MouseBrain_008um"),
                            header = TRUE)
doublet_info <- read.table(paste0("visium_hd_brain/Visium_HD_MouseBrain_008um/proportions_rctd_Visium_HD_MouseBrain_008um_doublet_info.tsv"),
                           header = TRUE)

classes_colors <- c("singlet"="forestgreen", "doublet_certain"="navyblue",
                    "doublet_uncertain"="orange", "reject"="red")

# Add row names to deconv_props
rownames(doublet_props) <- doublet_info$spot

# Read scrublet results
scrublet_res <- read.csv("visium_hd_brain/Visium_HD_MouseBrain_008um/scrublet_Visium_HD_MouseBrain_008um.csv",
                         header = TRUE) %>% rename(spot_id = X)

scrublet_res_subset <- scrublet_res %>% filter(spot_id %in% doublet_info$spot) %>% 
  mutate(predicted_doublet = ifelse(predicted_doublet == "True", "doublet_certain", "singlet"))

# Confusion matrix of scrublet and rctd results
conf_matrix <- caret::confusionMatrix(scrublet_res_subset$predicted_doublet %>% factor(levels=names(classes_colors)),
                                      doublet_info$spot_class %>% factor(levels=names(classes_colors)))


pheatmap::pheatmap(conf_matrix$table[1:2,] %>% `rownames<-`(c("singlet", "doublet")),
                   scale="row",
                 color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
                 cluster_rows = F, cluster_cols = FALSE,
                 main = "Scrublet (rows) vs RCTD (cols) predictions",
                 display_numbers=conf_matrix$table[1:2,],
                 fontsize_number=12, number_color="white",
                 angle_col = 0,
                 filename = "visium_hd_brain/plots/conf_matrix_scrublet_vs_rctd.png",
                 width=6, height=3)
############################
# Read in segmentation
intersection <- read.csv("visium_hd_brain/intersections.csv")
# Columns:
# - area = area of entire segmented nucleus
# - nucleus_area = area of nucleus intersected with the spot

# Plot histogram of nucleus area
ggplot(intersection %>% distinct(id, area), aes(x=area)) +
  geom_histogram(bins=50) +
  #geom_vline(xintercept=200, color="red", linetype="dashed")+
  theme_minimal()

# Boxplot
ggplot(intersection %>% distinct(id, area),
       aes(x=1, y=area)) +
  geom_boxplot() +
  #geom_hline(yintercept=101.66, color="red", linetype="dashed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

filter_nucleus_area <- FALSE
filter_val <- quantile(intersection %>% distinct(id, area) %>% pull(area), probs = c(0.25))
if (filter_nucleus_area){
  intersection <- intersection %>% filter(area > filter_val)
}

grid_area <- max(intersection$nucleus_area) # 853.56

# Bin nucleus area in increments of 25
intersection <- intersection %>% mutate(nucleus_area_bin =
                    cut(nucleus_area, breaks=c(seq(0, max(intersection$nucleus_area), by=25), Inf)))

# Check -> does more nuclei area correspond to the counts?
ggplot(intersection, aes(x=nucleus_area_bin, y=total_counts)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("visium_hd_brain/plots/boxplot_counts_by_nucleas_area.png",
       width=10, height=5, dpi=300, bg="white")

# Per grid_id, label it as single, double or multiple depending on the number of ids
intersection <- intersection %>% group_by(grid_id) %>%
  mutate(nucleus_count = n()) %>% ungroup()

intersection %>% distinct(grid_id, nucleus_count) %>%
  filter(nucleus_count <=2) %>% pull(nucleus_count) %>% table %>% prop.table

# Barplot - number of nuclei per bin
ggplot(intersection %>% filter(area > filter_val) %>% 
         distinct(grid_id, nucleus_count),
       aes(x=factor(nucleus_count))) +
  geom_bar(width=0.6, fill="black") +
  labs(title="Number of nuclei per bin", x="Number of nuclei") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(size = 11),
        panel.grid.major.x = element_blank())
ggsave("visium_hd_brain/plots/barplot_nucleus_count_per_bin.png",
       width=7, height=5, dpi=300, bg="white")

# Is there the same trend for more nuclei? -> yes
ggplot(intersection, aes(x=nucleus_area_bin, y=total_counts)) +
  geom_boxplot() +
  theme_minimal() +
  facet_wrap(~nucleus_count) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

nucleus_pct_threshold <- 0

intersection_join_rctd <- inner_join(intersection, doublet_info, by = c("spot_id" = "spot")) #%>%
  # Only keep if nucleus_area_percentage > threshold
  filter(nucleus_area_percentage > nucleus_pct_threshold) %>% 
  # Only keep grid_id where the remaining spots matches with nucleus_count
  group_by(grid_id) %>% filter(n() == unique(nucleus_count))

intersection_join_rctd <- intersection_join_rctd %>%
  distinct(grid_id, spot_id, spot_class, nucleus_count) 

# Bar plot of n_nucleus and count spot_class
barplot_abs <- ggplot(intersection_join_rctd %>% filter(nucleus_count <= 2),
       aes(x=factor(nucleus_count), fill=factor(spot_class, levels = names(classes_colors)))) +
  geom_bar(width=0.6) +
  scale_fill_manual(values=classes_colors, name="Spot class") +
  labs(title="Absolute number of spots", x="Number of nuclei") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(size = 11),
        panel.grid.major.x = element_blank()
        )

barplot_rel <- ggplot(intersection_join_rctd %>% filter(nucleus_count <= 2),
       aes(x=factor(nucleus_count), fill=factor(spot_class, levels=names(classes_colors)))) +
  geom_bar(position="fill", width=0.6) +
  scale_fill_manual(values=classes_colors, name="Spot class") +
  labs(title="Relative fraction of spots", x="Number of nuclei") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(size = 11),
        panel.grid.major.x = element_blank())

barplot_abs + barplot_rel +
  # gather legends
  plot_layout(guides = "collect")
ggsave("visium_hd_brain/plots/barplot_rctd_spotclass_per_n_nucleus.png",
       width=7, height=5, dpi=300, bg="white")

intersection_join_scrub <- inner_join(intersection, scrublet_res, by = c("spot_id" = "spot_id")) %>%
  # Only keep if nucleus_area_percentage > threshold
  filter(nucleus_area_percentage > nucleus_pct_threshold) %>% 
  # Only keep grid_id where the remaining spots matches with nucleus_count
  group_by(grid_id) %>% filter(n() == unique(nucleus_count)) 
  
intersection_join_scrub_summ <- intersection_join_scrub %>%
  distinct(grid_id, spot_id, predicted_doublet, nucleus_count)

# Bar plot of nucleus_count and count spot_class
barplot_abs <- ggplot(intersection_join_scrub_summ %>% filter(nucleus_count <= 2),
       aes(x=factor(nucleus_count), fill=factor(predicted_doublet, levels = c("False", "True"))) ) +
  geom_bar(width=0.6) +
  scale_fill_manual(values=c("False"="forestgreen", "True"="navyblue"),
                    labels=c("False"="Singlet", "True"="Doublet")) +
  labs(title="Absolute number of spots", x="Number of nuclei") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(size = 11),
        panel.grid.major.x = element_blank(),
        legend.title = element_blank()
        )
barplot_abs
barplot_rel <- ggplot(intersection_join_scrub_summ %>% filter(nucleus_count <= 2),
       aes(x=factor(nucleus_count), fill=factor(predicted_doublet, levels = c("False", "True"))) ) +
  geom_bar(position="fill", width=0.6) +
  scale_fill_manual(values=c("False"="forestgreen", "True"="navyblue"),
                    labels=c("False"="Singlet", "True"="Doublet")) +
  labs(title="Relative fraction of spots", x="Number of nuclei") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(size = 11),
        panel.grid.major.x = element_blank(),
        legend.title = element_blank())
barplot_rel

barplot_abs + barplot_rel +
  # gather legends
  plot_layout(guides = "collect")
ggsave("visium_hd_brain/plots/barplot_scrublet_per_n_nucleus.png",
       width=7, height=5, dpi=300, bg="white")

# Calculating AUPR
intersection_join_scrub_filter <- intersection_join_scrub %>% filter(nucleus_count <= 2)
curves <- evalmod(mmdata(intersection_join_scrub_filter %>% pull(doublet_score),
                         intersection_join_scrub_filter %>% pull(nucleus_count)))
# Get auc
print(paste0("AUROC: ", subset(auc(curves), curvetypes=="ROC")$aucs)) #0.521
print(paste0("AUPRC: ", subset(auc(curves), curvetypes=="PRC")$aucs)) #0.363

classes <- c("True", "False")
confusionMatrix(intersection_join_scrub_filter$predicted_doublet %>% factor(levels=classes),
                intersection_join_scrub_filter %>% mutate(ground_truth = ifelse(nucleus_count == 1, "False", "True")) %>% 
                pull(ground_truth) %>% factor(levels=classes),
                mode = "prec_recall")
