source("visium_hd_liver_combined/0_utils.R")

deconv_props_rank <- readRDS(paste0("visium_hd_liver_combined/rds/deconv_props_all_doubletmode.rds"))

dataset_name <- "sca002" # sca002 or caw009
dataset <- ifelse(dataset_name == "caw009", "_caw009", "")
data_path <- paste0("data/Visium_HD_Liver", toupper(dataset), "/")
proportions_path <- paste0("visium_hd_liver", dataset, "/Visium_HD_Liver", toupper(dataset))
bin_size <- 16
bin_size_str <- sprintf("%03dum", bin_size)
visium_obj <- readRDS(paste0(data_path, "Visium_HD_Liver", toupper(dataset), "_",
                             bin_size_str, ".rds"))

# We want to check correlation between markers and deconvolution results
marker_genes <- read.csv("visium_hd_liver/MERSCOPE_marker_genes_onehot.csv")
lfc <- read.csv("data/liver_marker_genes/final_MERscope_panel_Guilliams1.csv")
KC_markers <- marker_genes %>% 
  # row where KC column is 1
  filter(.[["KC"]] == 1) %>% select(Vizgen.Gene) %>% 
  inner_join(lfc %>% select(Vizgen.Gene, Abundance)) %>% 
  slice_max(Abundance, n=5) %>% pull(Vizgen.Gene)

visium_obj_subset <- visium_obj[KC_markers, ]

expr_df <- as.matrix(GetAssayData(visium_obj_subset, layer = "counts")) %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("spot") %>% 
  pivot_longer(-spot, names_to = "gene", values_to = "expression")
  
deconv_props_KC <- deconv_props_rank %>% filter(dataset==dataset_name, bin_size==bin_size_str, celltype == "Kupffercells")
deconv_props_KC %>% dim
sum(!Cells(visium_obj_subset) %in% unique(deconv_props_KC$spot))

expr_df_KC <- expr_df %>% 
  left_join(deconv_props_KC) %>% 
  # fill NAs with 0
  mutate(proportion = ifelse(is.na(proportion), 0, proportion))

# Calculate correlation between Clec4f and proportion
expr_df_KC %>% filter(gene == "Slc40a1") %>% 
  summarise(correlation = cor(expression, proportion, use = "complete.obs"))

expr_df_KC <-  expr_df_KC %>% 
  # bin expression into 10 bins, but make one bin for just 0
  mutate(expression_bin = cut(expression, 
                              breaks = c(-Inf, seq(0, max(expression, na.rm = TRUE), length.out = 10)), 
                              include.lowest = TRUE))

expr_df_KC$expression_bin

# Plot genes
ggplot(expr_df_KC) +
  geom_boxplot(aes(x = expression_bin, y = proportion)) +
  facet_wrap(~gene) +
  theme_minimal() +
  labs(title = paste("Expression vs Proportion for", dataset_name, "bin size", bin_size_str),
       x = "Gene expression",
       y = "Proportion")

ggplot(expr_df_KC) +
  geom_point(aes(x = expression, y = proportion)) +
  facet_wrap(~gene) +
  theme_minimal() +
  labs(title = paste("Expression vs Proportion for", dataset_name, "bin size", bin_size_str),
       x = "Gene expression",
       y = "Proportion")
