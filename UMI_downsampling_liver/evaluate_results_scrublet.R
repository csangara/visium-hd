library(tidyverse)
library(Matrix)
library(caret)
library(shadowtext)
library(precrec)

repl <- 2
# Analyze downsampling results
doublet_info <- read.csv(paste0("UMI_downsampling/brain_cortex_generation_real/scrublet_brain_cortex_generation_real_rep",
                                  repl, ".csv"), header = TRUE) %>% 
  rename(spot=X) %>% 
  filter(!is.na(doublet_score))

synth_obj <- readRDS(paste0("data/MouseBrain_UMI_downsampling/brain_cortex_generation_real_rep", repl, ".rds"))

# Get ground truth singlet and doublet
ground_truth_class <- synth_obj$spot_composition %>% 
  filter(name %in% doublet_info$spot) %>%
  mutate(n_cells = rowSums(across(Astro:VLMC)),
        spot_class = case_when(n_cells == 1 ~ "FALSE",
                               n_cells == 2 ~ "TRUE")) %>% 
  select(-(Astro:VLMC))

# Find current doublet_score threshold (~0.18)
doublet_info %>% filter(predicted_doublet == "True") %>% pull(doublet_score) %>% min
doublet_info %>% filter(predicted_doublet == "False") %>% pull(doublet_score) %>% max

ggplot(doublet_info, aes(x=doublet_score)) +
  geom_histogram() +
  geom_vline(xintercept = 0.18, color="red", linetype="dashed") +
  theme_classic()
ggsave("UMI_downsampling/plots/histogram_scrublet_doublet_score.png", width=8, height=5)

# Find a doublet_score threshold that predicts 70% of data as singlet
doublet_info %>% 
  mutate(predicted_doublet = ifelse(doublet_score > 0.7, "TRUE", "FALSE")) %>% 
  pull(predicted_doublet) %>% table %>% 
  prop.table

# Overwrite predicted_doublet column using 0.5 as threshold
doublet_info <- doublet_info %>% 
  mutate(predicted_doublet = ifelse(doublet_score > 0.7, "TRUE", "FALSE"))

classes <- c("TRUE", "FALSE")
confusionMatrix(doublet_info$predicted_doublet %>% factor(levels=classes),
                ground_truth_class$spot_class %>% factor(levels=classes),
                mode = "prec_recall")

# AUPR?
evalmod(mmdata(doublet_info$doublet_score,
       ground_truth_class$n_cells))

# Divide per bin of 100 UMIs
ground_truth_counts <- synth_obj$counts
gt_counts_df <- ground_truth_counts %>% colSums %>% data.frame(nCounts = .) %>% 
  mutate(counts_bin = cut(nCounts, breaks=c(seq(0, 2000, by=100), Inf))) %>%
  rownames_to_column("spot") %>% 
  inner_join(ground_truth_class, by=c("spot"="name"))

bin_names <- str_replace_all(levels(gt_counts_df$counts_bin), "e\\+03", "k") %>% 
  setNames(levels(gt_counts_df$counts_bin))
bin_names_right <- c(paste0("\u2264", seq(100, 2000, by=100)), ">2000") %>%
  setNames(levels(gt_counts_df$counts_bin))

# Calculate confusion matrix of singlet/doublet prediction
conf_matrices <- gt_counts_df %>%
  group_by(counts_bin) %>%
  group_split() %>%
  lapply(
    function(gt_count_df){
      confusionMatrix(doublet_info %>%
                        filter(spot %in% gt_count_df$spot) %>%
                        pull(predicted_doublet) %>% factor(levels=classes),
                      gt_count_df$spot_class %>% factor(levels=classes),
                      mode = "prec_recall")
    })


# Get accuracy
acc_df <- purrr::map(conf_matrices, "overall") %>%
  setNames(gt_counts_df$counts_bin %>% levels) %>% 
  do.call(rbind, .) %>% data.frame %>% 
  rownames_to_column("bins") %>% 
  mutate(bins = factor(bins, levels=names(bin_names))) %>% 
  inner_join(gt_counts_df %>% group_by(counts_bin) %>% summarise(n=n()),
             by =c("bins" = "counts_bin"))

# Get other classification metrics
class_metrics_df <- lapply(conf_matrices, function(mat) {
  mat$byClass %>% data.frame(value=.) %>%
    rownames_to_column("metric")
}) %>% setNames(names(bin_names)) %>% 
  list_rbind(., names_to = "bins") %>% 
  mutate(bins = factor(bins, levels=names(bin_names)))

ggplot(acc_df, aes(x=bins, y=Accuracy, group=1)) +
  geom_line() +
  geom_point(aes(size=n)) +
  scale_x_discrete(breaks=names(bin_names_right)[c(1, 5, 10, 15, 20, 21)],
                   labels=bin_names_right[c(1, 5, 10, 15, 20, 21)]) +
  theme_classic() +
  ylim(0, 1) +
  labs(x = "UMI Counts Per Spot", y="Accuracy of spot class prediction", size="# spots in bin") +
  theme(legend.position.inside=c(0.15, 0.9),
        legend.position="inside")
ggsave("UMI_downsampling/plots/scrublet_spotclass_accuracy_per_bin.png", width=8, height=8)

ggplot(class_metrics_df %>% filter(metric=="F1"),
       aes(x=bins, y=value, group=1)) +
  geom_line() +
  geom_point() +
  scale_x_discrete(breaks=names(bin_names_right)[c(1, 5, 10, 15, 20, 21)],
                   labels=bin_names_right[c(1, 5, 10, 15, 20, 21)]) +
  scale_color_discrete(labels=c("Singlet", "Doublet (certain)")) +
  theme_classic() +
  ylim(0, 1) +
  labs(x = "UMI Counts Per Spot", y="F1 Score", color="Predicted spot class") +
  theme(legend.position.inside=c(0.15, 0.9),
        legend.position="inside")
ggsave("UMI_downsampling/plots/scrublet_spotclass_F1_per_bin.png", width=8, height=8)


conf_matrices_df <-  purrr::map(conf_matrices, "table") %>% 
  lapply (function(mat) {
    mat %>% data.frame %>% 
      filter(!Reference %in% c("doublet_uncertain", "reject")) %>% 
      mutate(Prediction = factor(Prediction, levels=rev(classes)),
             Reference = factor(Reference, levels=classes),
             Freq = Freq/sum(Freq)*100)}) %>% 
  setNames(levels(gt_counts_df$counts_bin)) %>% 
  list_rbind(., names_to = "bins") %>% 
  mutate(bins = factor(bins, levels=levels(gt_counts_df$counts_bin)))


p_conf_mat <- ggplot(conf_matrices_df, aes(x=Reference, y=Prediction, fill=Freq)) +
  geom_tile() +
  #geom_shadowtext(data=conf_matrices_df %>% filter(round(Freq, 2) > 0),
  #                aes(label=paste0(round(Freq), "%"))) +
  facet_wrap(~bins, nrow = 1, labeller=labeller(bins = bin_names_right), strip.position = "bottom") +
  scale_fill_gradientn(colours=RColorBrewer::brewer.pal(9,"Blues")) +
  scale_x_discrete(position="top", labels=c("singlet", "doublet")) +
  labs(fill="Freq(%)") +
  theme_minimal() +
  coord_fixed() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, hjust=0),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_text(size=11))

p_conf_mat_tab <- ggplotGrob(p_conf_mat)

ggsave("UMI_downsampling/plots/scrublet_spotclass_confusion_matrix_per_bin.png", 
       # Filter out top axes 2-21 (only keep the first one)
       plot=grid::grid.draw(gtable::gtable_filter(p_conf_mat_tab, "axis-t-([2-9]|1[0-9]|2[0-1])", invert=TRUE)),
       width=20, height=3, bg="white")

# Are the percentages of spot class prediction different across bins?
doublet_info_bins <- doublet_info %>% 
  inner_join(gt_counts_df %>% select(spot, counts_bin), by="spot")

library(ggridges)
ggplot(doublet_info_bins, aes(x=doublet_score, y=forcats::fct_rev(counts_bin))) +
  scale_y_discrete(breaks=names(bin_names_right)[c(1, 5, 10, 15, 20, 21)],
                   labels=bin_names_right[c(1, 5, 10, 15, 20, 21)]) +
  geom_density_ridges() +
  labs(y="UMI Counts Per Spot", x="Doublet Score") +
  theme_classic()
ggsave("UMI_downsampling/plots/scrublet_doublet_score_per_bin.png", width=6, height=5)

# Stacked bar plot
ggplot(doublet_info_bins, aes(y=forcats::fct_rev(counts_bin), fill=factor(predicted_doublet, levels=classes))) +
  geom_bar(position="fill") +
  scale_fill_manual(values=c("FALSE"="forestgreen", "TRUE"="navyblue"),
                    labels=c("Singlet", "Doublet"),
                    breaks=c("FALSE", "TRUE")) +
  scale_y_discrete(breaks=names(bin_names_right)[c(1, 5, 10, 15, 20, 21)],
                    labels=bin_names_right[c(1, 5, 10, 15, 20, 21)]) +
  theme_classic() +
  labs(y = "UMI Counts Per Spot", x="Proportion of spots", fill="Predicted spot class")
ggsave("UMI_downsampling/plots/scrublet_spotclass_proportion_per_bin.png", width=8, height=5)

# Get AUPR per bin
gt_binary_bins <- gt_counts_df %>% select(spot, counts_bin, n_cells) %>% 
  column_to_rownames("spot") %>% 
  group_by(counts_bin) %>% 
  group_split(.keep=FALSE) %>% 
  purrr::map(pull, n_cells)

doublet_props_bins <- gt_counts_df %>% select(spot, counts_bin) %>% 
  inner_join(doublet_info %>% select(spot, doublet_score), by="spot") %>% 
  column_to_rownames("spot") %>% 
  group_by(counts_bin) %>% 
  group_split(.keep=FALSE) %>% 
  purrr::map(pull, doublet_score)

# Calculate doublet classification AUPR per bin
model <- mmdata(doublet_props_bins, gt_binary_bins, dsids=1:21)
curve <- evalmod(model)
prc <- subset(auc(curve), curvetypes=="PRC") %>% rename(method=modnames)
ggplot(prc, aes(x=factor(dsids), y=aucs)) +
  geom_point() +
  geom_line(aes(group=1)) +
  scale_x_discrete(breaks=c(1, 5, 10, 15, 20, 21),
                   labels=bin_names_right %>% setNames(1:21) %>%
                     .[c(1, 5, 10, 15, 20, 21)]) +
  ylim(c(0,1)) +
  labs(x = "UMI Counts Per Spot", y="Doublet classification AUPR") +
  theme_classic()

ggsave("UMI_downsampling/plots/scrublet_doublet_AUPR_per_bin.png", width=8, height=6)
