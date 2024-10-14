library(tidyverse)
library(Matrix)
library(caret)
library(shadowtext)
library(precrec)

repl <- 1
# Analyze downsampling results
doublet_props <- read.table(paste0("UMI_downsampling/brain_cortex_generation_real/proportions_nnls_brain_cortex_generation_real_rep",
                                   repl), header = TRUE)

synth_obj <- readRDS(paste0("data/MouseBrain_UMI_downsampling/brain_cortex_generation_real_rep", repl, ".rds"))
ground_truth <- synth_obj$relative_spot_composition

ground_truth %>% dim # 11266
ground_truth <- ground_truth %>% select(-region) %>% 
  column_to_rownames("name") %>%
  `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", "")) %>% 
  .[,sort(colnames(.), method="shell")]

# Correlation and RMSE
corr_spots <- mean(diag(cor(t(ground_truth), t(doublet_props))))
RMSE <- mean(sqrt(rowSums((ground_truth-doublet_props)**2)/ncol(ground_truth)))

# Load counts
ground_truth_counts <- synth_obj$counts

# Divide per bin of 100 UMIs
gt_counts_df <- ground_truth_counts %>% colSums %>% data.frame(nCounts = .) %>% 
  mutate(counts_bin = cut(nCounts, breaks=c(seq(0, 2000, by=100), Inf))) %>%
  rownames_to_column("spot")

bin_names <- str_replace_all(levels(gt_counts_df$counts_bin), "e\\+03", "k") %>% 
  setNames(levels(gt_counts_df$counts_bin))
bin_names_right <- c(paste0("\u2264", seq(100, 2000, by=100)), ">2000") %>%
  setNames(levels(gt_counts_df$counts_bin))

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
               mutate(spot = colnames(ground_truth_counts)), by="spot") %>% 
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


# Calculate AUPR
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
  labs(x = "UMI Counts Per Spot", y="Cell type prediction AUPR") +
  theme_classic()
ggsave("UMI_downsampling/plots/celltype_AUPR_per_bin.png", width=8, height=8)


