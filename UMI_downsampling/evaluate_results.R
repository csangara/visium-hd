library(tidyverse)
library(Matrix)
library(caret)
library(shadowtext)
library(precrec)
library(ggridges)

repl <- 2
# Analyze downsampling results
doublet_props <- read.table(paste0("UMI_downsampling/brain_cortex_generation_real/proportions_rctd_brain_cortex_generation_real_rep",
                                   repl, "_2.2.1-hd-seurat5.1.0"), header = TRUE)
doublet_info <- read.table(paste0("UMI_downsampling/brain_cortex_generation_real/proportions_rctd_brain_cortex_generation_real_rep",
                                  repl, "_2.2.1-hd-seurat5.1.0_doublet_info.tsv"), header = TRUE)

synth_obj <- readRDS(paste0("data/MouseBrain_UMI_downsampling/brain_cortex_generation_real_rep", repl, ".rds"))
ground_truth <- synth_obj$relative_spot_composition

# Subset ground truth to rows in doublet_info
ground_truth %>% dim # 11266
ground_truth <- ground_truth %>% filter(name %in% doublet_info$spot) %>% select(-region) %>% 
  column_to_rownames("name") %>%
  `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", "")) %>% 
  .[,sort(colnames(.), method="shell")]

# Correlation and RMSE
corr_spots <- mean(diag(cor(t(ground_truth), t(doublet_props))))
RMSE <- mean(sqrt(rowSums((ground_truth-doublet_props)**2)/ncol(ground_truth)))

# Get ground truth singlet and doublet
ground_truth_class <- synth_obj$spot_composition %>% 
  filter(name %in% doublet_info$spot) %>% 
  mutate(n_cells = rowSums(across(Astro:VLMC)),
        spot_class = case_when(n_cells == 1 ~ "singlet",
                               n_cells == 2 ~ "doublet_certain")) %>% 
  select(-(Astro:VLMC))

classes <- doublet_info$spot_class %>% unique

confusionMatrix(doublet_info$spot_class %>% factor(levels=classes),
                ground_truth_class$spot_class %>% factor(levels=classes))

# Load counts
ground_truth_counts <- synth_obj$counts %>%
  .[, colnames(.) %in% doublet_info$spot]

# Histogram
png(paste0("UMI_downsampling/plots/synth_obj_rep", repl,"_counts_hist.png"), width=6, height=4, units="in", res=300)
ground_truth_counts %>% colSums %>% hist(breaks=100, xlab="UMI Counts", main=paste0("Histogram of rep ", repl))
dev.off()

# Compare to silver standard
# png(paste0("UMI_downsampling/plots/silver_standard_1-10_rep1_counts_hist_log.png"), width=6, height=4, units="in", res=300)
# readRDS("~/spotless-benchmark/standards/silver_standard_1-10/brain_cortex_real_rep1.rds")$counts %>% colSums %>% 
#   hist(breaks=100, xlab="UMI Counts", main="Histogram of silver standard")
# dev.off()


# Divide per bin of 100 UMIs
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
                    pull(spot_class) %>% factor(levels=classes),
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
  mat$byClass %>% data.frame() %>%
    rownames_to_column("type") %>%
    mutate(type = str_remove(pattern = "Class: ", type)) %>%
    pivot_longer(-type)
}) %>% setNames(names(bin_names)) %>% 
  list_rbind(., names_to = "bins") %>% 
  mutate(bins = factor(bins, levels=names(bin_names)))

ggplot(acc_df, aes(x=bins, y=Accuracy, group=1)) +
  geom_line() +
  geom_point(aes(size=n)) +
  scale_x_discrete(breaks=names(bin_names_right)[c(1, 5, 10, 15, 20, 21)],
                   labels=bin_names_right[c(1, 5, 10, 15, 20, 21)]) +
  theme_classic(base_size=24) +
  ylim(0, 1) +
  labs(x = "UMI Counts Per Spot", y="Accuracy of spot class prediction", size="# spots in bin") +
  theme(legend.position.inside=c(0.25, 0.9),
        legend.position="inside")
ggsave("UMI_downsampling/plots/spotclass_accuracy_per_bin.png", width=8, height=8, dpi=300)

ggplot(class_metrics_df %>% filter(name=="F1", type %in% c("singlet", "doublet_certain")),
       aes(x=bins, y=value, color=factor(type, levels=c("singlet", "doublet_certain")), group=type)) +
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
ggsave("UMI_downsampling/plots/spotclass_F1_per_bin.png", width=8, height=8)

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


ggsave("UMI_downsampling/plots/spotclass_confusion_matrix_per_bin.png", 
       # Filter out top axes 2-21 (only keep the first one)
       plot=grid::grid.draw(gtable::gtable_filter(p_conf_mat_tab, "axis-t-([2-9]|1[0-9]|2[0-1])", invert=TRUE)),
       width=20, height=3, bg="white")

# Can min score be used to predict doublets? -> no
ggplot(doublet_info, aes(x=min_score, y=spot_class)) +
  geom_density_ridges() +
  theme_classic()

ggplot(doublet_info, aes(x=singlet_score, y=min_score, color=factor(spot_class, levels=classes))) +
  geom_point(size=0.01) +
  coord_fixed() +
  scale_color_manual(values=c("singlet"="forestgreen", "doublet_certain"="navyblue",
                             "doublet_uncertain"="orange", "reject"="red"),
                     name="Doublet class") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  #facet_wrap(~spot_class) +
  theme_classic()
ggsave("UMI_downsampling/plots/singlet_vs_min_score.png", width=8, height=5)

gt_binary_bins <- gt_counts_df %>% select(spot, counts_bin, n_cells) %>%
  filter(spot %in% doublet_info$spot) %>%
  column_to_rownames("spot") %>% 
  group_by(counts_bin) %>% 
  group_split(.keep=FALSE) %>% 
  purrr::map(pull, n_cells)

doublet_props_bins <- gt_counts_df %>% select(spot, counts_bin) %>% 
  inner_join(doublet_info %>% select(spot, min_score), by="spot") %>% 
  column_to_rownames("spot") %>% 
  group_by(counts_bin) %>% 
  group_split(.keep=FALSE) %>% 
  purrr::map(pull, min_score)

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

# Are the percentages of spot class prediction different across bins?
doublet_info_bins <- doublet_info %>% 
  inner_join(gt_counts_df %>% select(spot, counts_bin), by="spot")

# Stacked bar plot
ggplot(doublet_info_bins, aes(y=forcats::fct_rev(counts_bin), fill=factor(spot_class, levels=rev(classes)))) +
  geom_bar(position="fill") +
  scale_fill_manual(values=c("singlet"="forestgreen", "doublet_certain"="navyblue",
                             "doublet_uncertain"="orange", "reject"="red"),
                    breaks=classes) +
  scale_y_discrete(breaks=names(bin_names_right)[c(1, 5, 10, 15, 20, 21)],
                    labels=bin_names_right[c(1, 5, 10, 15, 20, 21)]) +
  theme_classic() +
  labs(y = "UMI Counts Per Spot", x="Proportion of spots", fill="Predicted spot class")
ggsave("UMI_downsampling/plots/spotclass_proportion_per_bin.png", width=8, height=5)

# Check AUPR per bin
gt_binary_bins <- gt_counts_df %>% select(spot, counts_bin) %>% 
  inner_join(ifelse(ground_truth > 0, "present", "absent") %>% data.frame %>% 
               rownames_to_column("spot"), by="spot") %>% 
  column_to_rownames("spot") %>% 
  group_by(counts_bin) %>% 
  group_split(.keep=FALSE) %>% 
  purrr::map(function(x) c(as.matrix(x)))

doublet_props_bins <- gt_counts_df %>% select(spot, counts_bin) %>% 
  inner_join(doublet_props %>% data.frame %>% 
               mutate(spot = doublet_info$spot), by="spot") %>% 
  column_to_rownames("spot") %>% 
  group_by(counts_bin) %>% 
  group_split(.keep=FALSE) %>% 
  purrr::map(function(x) c(as.matrix(x)))

celltypes <- colnames(doublet_props)

# Calculate overall AUPR
gt_binary <- ifelse(ground_truth > 0, "present", "absent") %>%
  reshape2::melt() %>% dplyr::select(value)

# Area under precision-recall curve
eval_prc <- evalmod(scores = c(as.matrix(doublet_props)), labels=gt_binary)
prc <- subset(auc(eval_prc), curvetypes == "PRC")$aucs
prc

# Calculate AUPR per bin
model <- mmdata(doublet_props_bins, gt_binary_bins, dsids=1:21)
curve <- evalmod(model)
prc <- subset(auc(curve), curvetypes=="PRC") %>% rename(method=modnames)
ggplot(prc, aes(x=factor(dsids), y=aucs)) +
  geom_point(size=3) +
  geom_line(aes(group=1)) +
  scale_x_discrete(breaks=c(1, 5, 10, 15, 20, 21),
                   labels=bin_names_right %>% setNames(1:21) %>%
                     .[c(1, 5, 10, 15, 20, 21)]) +
  ylim(c(0,1)) +
  labs(x = "UMI Counts Per Spot", y="Cell type prediction AUPR") +
  theme_classic(base_size=24)
ggsave("UMI_downsampling/plots/celltype_AUPR_per_bin.png", width=8, height=8)

# What is the AUPR of cell type prediction per spot class per bin?
aupr_curves <- lapply(levels(gt_counts_df$counts_bin), function(bin) {
  props_list <- lapply(c("ground_truth", "predicted"), function(x) {
    gt_counts_df %>% 
      filter(counts_bin %in% bin) %>% 
      select(spot) %>% 
      # Get spot prediction
      inner_join(doublet_info %>% select(spot, spot_class), by="spot") %>% 
      # Join with corresponding proportions (ground truth or predicted)
      {if (x == "ground_truth") {
        inner_join(., ifelse(ground_truth > 0, "present", "absent") %>% data.frame %>% 
                     rownames_to_column("spot"), by="spot")
      } else {
        inner_join(., doublet_props %>% data.frame %>% 
                     mutate(spot = doublet_info$spot), by="spot")
      }} %>%
      column_to_rownames("spot") %>% 
      mutate(spot_class = factor(spot_class, levels=classes)) %>% 
      group_by(spot_class) %>% 
      group_split(.keep=FALSE) %>% 
      purrr::map(function(x) c(as.matrix(x)))
  
  })

  model <- mmdata(props_list[[2]], props_list[[1]], dsids=1:4)
  curve <- evalmod(model)
  prc <- subset(auc(curve), curvetypes=="PRC") %>% rename(method=modnames) %>% 
    select(-curvetypes) %>% mutate(bin=bin)
  prc
}) %>% do.call(rbind, .)

ggplot(aupr_curves, aes(x=factor(dsids), y=aucs)) +
  geom_boxplot() +
  ylim(c(0,1)) +
  scale_x_discrete(labels=classes) +
  ggtitle("AUPR of cell type classification") +
  theme_classic() +
  theme(axis.title = element_blank())
ggsave("UMI_downsampling/plots/boxplot_celltype_AUPR_per_spotclass_per_bin.png", width=6, height=4)

