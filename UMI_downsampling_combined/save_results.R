library(tidyverse)
library(Matrix)
library(caret)
library(precrec)
library(patchwork)

get_file_paths <- function(repl, dataset, tool){
  directory_ext <- ifelse(dataset == "liver", "_liver", "_brain")
  file_prefix <- ifelse(dataset == "liver", "liver_mouseStSt_exVivo_9celltypes_annot_cd45_artificial_dominant_celltype_diverse", "brain_cortex_generation_real")
  file_suffix <- ifelse(dataset == "liver", "", "_2.2.1-hd-seurat5.1.0")
  
  if (tool == "rctd"){
    doublet_props_file <- paste0("UMI_downsampling", directory_ext, "/", file_prefix, "/proportions_rctd_", file_prefix, "_rep", repl, file_suffix)
    doublet_info_file <- paste0(doublet_props_file, "_doublet_info.tsv")
  } else if (tool == "scrublet") {
    doublet_props_file <- NA
    doublet_info_file <- paste0("UMI_downsampling", directory_ext, "/", file_prefix, "/scrublet_", file_prefix, "_rep", repl, ".csv")
  }
  synth_object <- paste0("UMI_downsampling_", dataset, "/", "synthetic_data/", file_prefix, "_rep", repl, ".rds")
  
  return(c(doublet_props_file,
           doublet_info_file,
           synth_object))
}


all_metrics <-
lapply(c("liver", "brain"), function(ds){
lapply(c("rctd", "scrublet"), function(tool) {
  if (tool == "scrublet" & ds == "liver") {
    return (data.frame())
  }
lapply(1:3, function(k) {
  file_paths <- get_file_paths(repl = k, dataset = ds, tool = tool)
  
  if (tool == "rctd") {
    doublet_props <- read.table(file_paths[1], header = TRUE)
  } else {
    doublet_props <- NA
  }
  doublet_info <-  read.table(file_paths[2], header = TRUE,
                              sep = ifelse(tool == "scrublet", ",", "\t"))
  
  if (tool == "scrublet"){
    doublet_info <- doublet_info %>% rename(spot = X) %>% 
      filter(!is.na(doublet_score))
    
    # Find a doublet_score threshold that will return half/half predictions
    threshold <- quantile(doublet_info$doublet_score, 0.5)
    
    # doublet_info <- doublet_info %>% 
    #   mutate(spot_class = ifelse(predicted_doublet == "True", "doublet", "singlet"))
    
    # Overwrite predicted_doublet column using new htreshold
    doublet_info <- doublet_info %>%
      mutate(spot_class = ifelse(doublet_score > threshold, "doublet_certain", "singlet"))
    
  }

  synth_obj <- readRDS(file_paths[3])
  
  ground_truth <- synth_obj$relative_spot_composition
  
  # Subset ground truth to rows in doublet_info
  ground_truth <- ground_truth %>% filter(name %in% doublet_info$spot) %>%
    select(-region) %>% 
    column_to_rownames("name") %>%
    `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", "")) %>% 
    .[,sort(colnames(.), method="shell")]
  
  # Load counts
  ground_truth_counts <- synth_obj$counts %>%
    .[, colnames(.) %in% doublet_info$spot]
  
  # Get ground truth singlet and doublet
  ground_truth_class <- synth_obj$spot_composition %>% 
    filter(name %in% doublet_info$spot) %>% 
    mutate(n_cells = rowSums(across(!c(name, region))),
           spot_class = case_when(n_cells == 1 ~ "singlet",
                                  n_cells == 2 ~ "doublet_certain")) %>% 
    select(name:spot_class)
  
  # Divide per bin of 100 UMIs
  gt_counts_df <- ground_truth_counts %>% colSums %>% data.frame(nCounts = .) %>% 
    mutate(counts_bin = cut(nCounts, breaks=c(seq(0, 2000, by=100), Inf))) %>%
    rownames_to_column("spot") %>% 
    inner_join(ground_truth_class, by=c("spot"="name"))
  
  bin_names <- str_replace_all(levels(gt_counts_df$counts_bin), "e\\+03", "k") %>% 
    setNames(levels(gt_counts_df$counts_bin))
  bin_names_right <- c(paste0("\u2264", seq(100, 2000, by=100)), ">2000") %>%
    setNames(levels(gt_counts_df$counts_bin))
  
  classes <- doublet_info$spot_class %>% unique
  
  # classes <- unique(str_extract(doublet_info$spot_class, "[a-z]+"))
  # Calculate confusion matrix of singlet/doublet prediction
  conf_matrices <- gt_counts_df %>%
    group_by(counts_bin) %>%
    group_split() %>%
    lapply(
      function(gt_count_df){
        confusionMatrix(doublet_info %>%
                          filter(spot %in% gt_count_df$spot) %>%
                          #mutate(spot_class = str_extract(spot_class, "[a-z]+")) %>% 
                          pull(spot_class)  %>% factor(levels=classes),
                        gt_count_df$spot_class %>% factor(levels=classes),
                        mode = "prec_recall")
      })
  
  acc_df <- purrr::map(conf_matrices, "overall") %>%
    setNames(gt_counts_df$counts_bin %>% levels) %>% 
    do.call(rbind, .) %>% data.frame %>% 
    rownames_to_column("bins") %>% 
    mutate(bins = factor(bins, levels=names(bin_names))) %>% 
    inner_join(gt_counts_df %>% group_by(counts_bin) %>% summarise(n=n()),
               by =c("bins" = "counts_bin"))
  
  conf_matrices_df <-  purrr::map(conf_matrices, "table") %>% 
    lapply (function(mat) {
      mat %>% data.frame %>% 
        filter(!Reference %in% c("doublet_uncertain", "reject")) %>% 
        mutate(Prediction = factor(Prediction, levels=rev(classes)),
               Reference = factor(Reference, levels=classes))}) %>% 
    setNames(levels(gt_counts_df$counts_bin)) %>% 
    list_rbind(., names_to = "bins") %>% 
    mutate(bins = factor(bins, levels=levels(gt_counts_df$counts_bin)))
  
  
  # Are the percentages of spot class prediction different across bins?
  doublet_info_bins <- doublet_info %>% 
    inner_join(gt_counts_df %>% select(spot, counts_bin), by="spot") %>% 
    count(spot_class, counts_bin)
  
  if (tool == "rctd"){
    # Cell type AUPR
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
      group_split(.keep=FALSE)
    
    doublet_props_bins_long <- doublet_props_bins %>% 
      purrr::map(function(x) c(as.matrix(x)))
    
    celltypes <- colnames(doublet_props)
    
    # Calculate AUPR per bin
    model <- mmdata(doublet_props_bins_long, gt_binary_bins, dsids=1:21)
    curve <- evalmod(model)
    prc <- subset(auc(curve), curvetypes=="PRC") %>% rename(method=modnames)
    
    # What is the AUPR of cell type prediction per spot class?
    prcs_spot_classes <- lapply(levels(gt_counts_df$counts_bin), function(bin) {
      props_list <- lapply(c("ground_truth", "predicted"), function(x) {
        tmp <- gt_counts_df %>% 
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
          column_to_rownames("spot")
        classes <- unique(tmp$spot_class)
        tmp <- tmp %>% 
          mutate(spot_class = factor(spot_class, levels=classes)) %>% 
          group_by(spot_class) %>% 
          group_split(.keep=FALSE) %>% 
          purrr::map(function(x) c(as.matrix(x))) %>% setNames(classes)
        
      })
      classes <- names(props_list[[1]])
      model <- mmdata(props_list[[2]], props_list[[1]], dsids=1:length(classes))
      curve <- evalmod(model)
      prc <- subset(auc(curve), curvetypes=="PRC") %>% rename(method=modnames) %>% 
        select(-curvetypes) %>% mutate(bin=bin)
      prc %>% mutate(spot_class = classes)
      
    }) %>% do.call(rbind, .)
    
    df <- bind_rows(prc %>%
                      select(-c(method, dsids, curvetypes)) %>% 
                       mutate(metric = "AUPR", bins = names(bin_names), eval = "celltype_classification") %>% 
                       rename(value = aucs),
                    acc_df %>% select(bins, Accuracy) %>% 
                       mutate(metric = "Accuracy", eval = "doublet_classification") %>% 
                       rename(value = Accuracy),
                    prcs_spot_classes %>%
                      select(-c(method, dsids)) %>% 
                      rename(value = aucs, bins=bin) %>% 
                      mutate(metric = "AUPR", eval = "AUPR_spotclass"))
  } else {
    
    df <-  acc_df %>% select(bins, Accuracy) %>% 
      rename(value = Accuracy) %>% 
      mutate(metric = "Accuracy", eval = "doublet_classification")
    
  }
  
  df <- bind_rows(df,
                  doublet_info_bins %>% mutate(metric = "counts") %>% 
                    rename(value = n, bins=counts_bin) %>%
                    mutate(eval = "doublet_frac_per_bin"),
                  conf_matrices_df %>% mutate(metric = "freq") %>% 
                    rename(value = Freq) %>% 
                    mutate(eval = "confusion_matrix")
                  )
                
  return(df %>% mutate(tool = tool, dataset = ds, repl = k))

}) %>% bind_rows()
}) %>% bind_rows()
}) %>% bind_rows()


saveRDS(all_metrics, "UMI_downsampling_combined/all_metrics.rds")
all_metrics <- readRDS("UMI_downsampling_combined/all_metrics.rds")
all_metrics <- all_metrics %>% 
  mutate(bins = factor(bins,  levels=unique(all_metrics$bins)))

bin_names <- str_replace_all(levels(all_metrics$bins), "e\\+03", "k") %>% 
  setNames(levels(all_metrics$bins))
bin_names_right <- c(paste0("\u2264", seq(100, 2000, by=100)), ">2000") %>%
  setNames(levels(all_metrics$bins))

spotclass_colors <- c("singlet"="#4DAF4A", "doublet_certain"="#377EB8",
                           "doublet_uncertain"="#FF7F00", "reject"="#E41A1C")
proper_spotclass_names <- c("Singlet", "Doublet\ncertain", "Doublet\nuncertain", "Reject") %>% 
  setNames(names(spotclass_colors))
tool_colors <- c("rctd" = "#1B9E77", "scrublet" = "#7570B3")
bin_names_sub <- c("(0,100]", "(500,\n600]", "(1000,\n1100]", "(1500,\n1600]", "(2000,\nInf]")

plot_a <- ggplot(all_metrics %>%
                   filter(dataset != "liver", eval == "celltype_classification"),
                 aes(x=bins, y = value, group=interaction(repl, tool))) +
  geom_line(color = "#1B9E77", linewidth=0.25) +
  geom_point(fill = "white", color = "#1B9E77", size=0.6, shape=21, stroke=0.25) +
  scale_x_discrete(breaks=names(bin_names[seq(1, 21, by = 5)]),
                   labels=bin_names_sub) +
  ylim(c(0.4,1)) +
  labs(x = "UMI Counts Per Spot", y="Cell type classification AUPR") +
  theme_classic(base_size=7) +
  ggtitle("Cell type classification AUPR\nacross three replicates") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size=6),
        panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        panel.grid.major.x = element_line(),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))
plot_a


plot_b <- ggplot(all_metrics %>%
         filter(dataset != "liver", eval == "AUPR_spotclass") %>% 
           mutate(spot_class = factor(spot_class, levels = names(spotclass_colors))),
       aes(x=spot_class, y = value, group = interaction(repl, spot_class), color=spot_class)) +
  scale_color_manual(values=spotclass_colors) +
  geom_boxplot(width=0.4, linewidth=0.25, outlier.size=0.3) +
  theme_classic(base_size=7) +
  scale_x_discrete(labels=proper_spotclass_names) +
  ggtitle("AUPR of cell type classification\nper predicted doublet type") +
  theme(axis.title = element_blank(),
        panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        legend.position = "none")
  
plot_b

plot_c <- ggplot(all_metrics %>%
                   filter(dataset != "liver", eval == "doublet_classification"),
                 aes(x=bins, y = value, group=interaction(repl, tool), color = tool)) +
  geom_line(linewidth=0.25) +
  geom_point(fill = "white", size=0.6, shape=21, stroke=0.25) +
  scale_x_discrete(breaks=names(bin_names[seq(1, 21, by = 5)]),
                   labels=bin_names_sub) +
  #ylim(c(0.4,1)) +
  labs(x = "UMI Counts Per Spot", y="Cell type classification AUPR") +
  theme_classic(base_size=7) +
  ggtitle("Doublet classification accuracy\nacross three replicates") +
  scale_color_manual(values=tool_colors, name = "Tool", labels = c("RCTD", "Scrublet")) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size=6),
        panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        panel.grid.major.x = element_line(),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.2),
        legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(color="black", linewidth = 0.2))
plot_c

plot_d <- ggplot(all_metrics %>% filter(dataset != "liver", eval == "doublet_frac_per_bin") %>% 
                   mutate(spot_class = factor(spot_class, levels = names(spotclass_colors))),
                aes(x = value, y=forcats::fct_rev(bins), fill=forcats::fct_rev(spot_class))) +
  geom_bar(position="fill", stat = "identity") +
  scale_fill_manual(values=spotclass_colors, breaks=names(spotclass_colors),
                    labels=c("Singlet", "Doublet (certain)", "Doublet (uncertain)", "Reject")) +
  scale_y_discrete(breaks=names(bin_names[seq(1, 21, by = 5)]),
                   labels=bin_names_sub) +
  scale_x_continuous(labels=c(0, "", 0.5, "", 1),
                     expand = expansion(mult = c(0, 0.05))) +
  theme_bw(base_size=7) +
  facet_wrap(~tool, labeller = labeller(tool = c("rctd" = "RCTD", "scrublet" = "Scrublet"))) +
  labs(y = "UMI Counts Per Spot", x="Proportion of spots", fill="Spot class") +
  theme(panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth=0.2),
        panel.grid = element_blank(),
        axis.title = element_text(size=6),
        strip.background = element_blank(),
        legend.key.size = unit(0.3, "cm"))

plot_d

conf_df <- all_metrics %>% filter(dataset != "liver", eval == "confusion_matrix",
                       tool == "rctd") %>% group_by(inter = interaction(Prediction, Reference)) %>% 
                                                    #low = bins %in% names(bin_names[1:10])) %>% 
  summarise(value = sum(value)) %>% 
  separate(inter, into = c("Prediction", "Reference"), sep="\\.") %>% 
  mutate(Prediction = factor(Prediction, levels = rev(names(spotclass_colors))),
         Reference = factor(Reference, levels = names(spotclass_colors)))


plot_conf <- ggplot(conf_df,
       aes(x=Reference, y=Prediction, fill=value)) +
  geom_tile() +
  # geom_shadowtext(data=conf_matrices_df %>% filter(round(Freq, 2) > 0),
  #                aes(label=paste0(round(Freq), "%"))) +
  # facet_wrap(~low, strip.position = "bottom") +
  scale_fill_gradientn(colours=RColorBrewer::brewer.pal(9,"Blues"),
                       breaks = c(5000, 10000, 15000), labels = c("5,000", "", "15,000")) +
  scale_x_discrete(position="top", labels=c("Singlet", "Doublet")) +
  scale_y_discrete(labels=rev(c("Singlet", "Doublet\n(certain)", "Doublet\n(uncertain)", "Reject"))) +
  labs(fill="Frequency") +
  theme_minimal(base_size=7) +
  coord_fixed() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        #axis.text.x = element_text(angle=45, hjust=0),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.title.position = "top")
plot_conf


first_row <- plot_a + plot_c + free(plot_conf) + plot_layout(widths=c(1, 1, 0.5))
second_row <- plot_b + plot_d
all_plots <- first_row / second_row + plot_layout(heights = c(1, 1.2)) &
  theme(plot.title=element_text(size=7, hjust=0.5),
        axis.text = element_text(size=5.5))

#all_plots
ggsave("UMI_downsampling_combined/plots/umi_downsampling_all_plots.pdf",
       all_plots,
       width=150, height=140, units="mm")



counts_features_df <- lapply(1:3, function(k) {
  file_paths <- get_file_paths(repl = k, dataset = "brain", tool = "rctd")
  synth_obj <- readRDS(file_paths[3])
  
  seurat_obj <- Seurat::CreateSeuratObject(synth_obj$counts)
  
  bind_rows(
    data.frame(x = seurat_obj$nCount_RNA, 
               repl = k, type = "Count"),
    data.frame(x = seurat_obj$nFeature_RNA, 
               repl = k, type = "Feature")
  )
}) %>% bind_rows()

set1_colors <- RColorBrewer::brewer.pal(9, "Set1")[1:9]
set1_colors

# Compare counts of the two datasets
p_counts_features <- ggplot(counts_features_df, aes(x = repl, y=x, group=repl)) +
  geom_violin(linewidth = 0.25, fill="#999999") +
  geom_boxplot(width=0.1, linewidth = 0.25, outliers = FALSE) +
  facet_wrap(~type, scales = "free",
             labeller = labeller( .multi_line=FALSE)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x="Replicate") +
  theme_classic(base_size = 7) +
  theme(legend.position = "bottom",
        axis.title.x = element_text(size=9),
        axis.text = element_text(size = 9),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 10))
p_counts_features
ggsave(paste0("UMI_downsampling_combined/plots/umi_downsampling_histogram_counts_features.pdf"),
       p_counts_features,
       width = 8, height = 4)
counts_features_df$x %>% min

