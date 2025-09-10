library(Seurat)
library(tidyverse)
library(patchwork)

plot_path <- paste0("visium_hd_liver_combined/plots/")
bin_sizes <- c(8, 16, 32)


deconv_props_rank <- readRDS(paste0("visium_hd_liver_combined/rds/deconv_props_all.rds"))

# for (dataset_name in c("sca002", "caw009")){
#   dataset <- ifelse(dataset_name == "caw009", "_caw009", "")
#   data_path <- paste0("data/Visium_HD_Liver", toupper(dataset), "/")
#   proportions_path <- paste0("visium_hd_liver", dataset, "/Visium_HD_Liver", toupper(dataset))
#   
#   for (bin_size in bin_sizes){
#     gc()
#     bin_size_str <- sprintf("%03dum", bin_size)
#     
#     visium_obj <- readRDS(paste0(data_path, "Visium_HD_Liver", toupper(dataset), "_",
#                                  bin_size_str, ".rds"))
#   }
# }

dataset_name <- "caw009"
dataset <- ifelse(dataset_name == "caw009", "_caw009", "")
data_path <- paste0("data/Visium_HD_Liver", toupper(dataset), "/")
proportions_path <- paste0("visium_hd_liver", dataset, "/Visium_HD_Liver", toupper(dataset))
bin_sizes <- c(16,32)
n_bins <- (bin_sizes[2]/bin_sizes[1])**2
bin_size_strs <- sprintf("%03dum", bin_sizes)

visium_obj <- readRDS(paste0(data_path, "Visium_HD_Liver", toupper(dataset), "_",
                             bin_size_strs[2], ".rds"))

counts_df <- visium_obj@meta.data[,paste0("nCount_Spatial.",bin_size_strs[2]),drop=FALSE] %>% 
  rownames_to_column("spot") %>%
  # Rename final column
  rename(nCount = paste0("nCount_Spatial.",bin_size_strs[2]))

# Read parquet file
tissue_positions <- arrow::read_parquet(
  paste0("visium_hd_liver_combined/rds/tissue_positions_", bin_size_strs[1], "_", bin_size_strs[2],
          toupper(dataset), ".parquet")) %>% 
  select(barcode, containing_barcode) %>% 
  # Filter to those with 4 barcodes per containing barcode
  group_by(containing_barcode) %>%
  filter(n() == n_bins) %>% ungroup()

deconv_props_subset <- deconv_props_rank %>%
  filter(bin_size %in% bin_size_strs, dataset == dataset_name, celltype != "Hepatocytes")

deconv_KCs <- deconv_props_subset %>%
  filter(celltype == "Kupffercells")

mother_bins <- tissue_positions %>% 
  select(barcode, containing_barcode) %>% 
  left_join(deconv_KCs %>% filter(bin_size == bin_size_strs[1]) %>%
              ungroup() %>% select(spot, proportion) %>%
               rename(prop_child = proportion),
             by = c("barcode" = "spot")) %>% 
  inner_join(deconv_KCs %>% filter(bin_size == bin_size_strs[2]) %>%
              ungroup() %>% select(spot, proportion) %>% 
              rename(prop_mother = proportion),
            by = c("containing_barcode" = "spot")) %>% 
  # replace prop by 0 if it's NA
  mutate(prop_child = ifelse(is.na(prop_child), 0, prop_child))

mother_bins_range <- mother_bins %>% 
  group_by(containing_barcode) %>% 
  mutate(range = max(prop_child) - min(prop_child)) %>% 
  inner_join(counts_df, by = c("containing_barcode" = "spot"))

# ggplot(mother_bins_range) +
#   geom_point(aes(x=range, y=prop_mother), alpha=0.5)

ggplot(mother_bins_range %>% ungroup() %>% 
         distinct(containing_barcode, prop_mother, range, nCount) %>% 
           arrange(nCount), 
       aes(x=prop_mother, y=range, color=nCount)) +
  geom_abline(slope=c(1), intercept=0, color="gray70", linetype="dashed") +
  #geom_abline(slope=2, intercept=0, color="gray70", linetype="dashed") +
  geom_point(size=0.5) +
  labs(x=paste0("Proportions of Kupffer cells in mother bins (", bin_sizes[2], "um)"),
       y=paste0("Range of proportions across ", n_bins, " child bins (", bin_sizes[1], "um)"),
       color="UMI counts\nin mother bin") +
  # Geom ab line
  scale_color_viridis_c() +
  theme_classic()
          #ggpointdensity::geom_pointdensity()

# Bin props_mother into 20 bins
mother_bins_binned <- mother_bins %>% 
  mutate(prop_mother_bin = cut(prop_mother, breaks=seq(0, 1, by=0.05), include.lowest=TRUE)) %>% 
  group_by(prop_mother_bin) %>% 
  mutate(bin_size = n()) %>% ungroup()

p <- ggplot(mother_bins_binned,
       aes(x=prop_mother_bin, y=prop_child)) +
  geom_boxplot()
ggsave(paste0(plot_path, "mother_bins_boxplot_kupffer_cells_", dataset_name, "_",
              bin_size_strs[1], "_", bin_size_strs[2], ".png"),
       p, width=8, height=4)
  
mean_prop <- deconv_KCs %>% filter(bin_size == bin_size_strs[2]) %>% 
  pull(proportion) %>% mean
var_prop <- deconv_KCs %>% filter(bin_size == bin_size_strs[2]) %>% 
  pull(proportion) %>% var

alpha <-mean_prop*((mean_prop*(1-mean_prop)/var_prop) - 1)
beta <- (1-mean_prop)*((mean_prop*(1-mean_prop)/var_prop) - 1)

sim_data <- rbeta(mother_bins %>% distinct(containing_barcode) %>% nrow(),
                    shape1 = alpha, shape2 = beta)

var_prop_child <- deconv_KCs %>% filter(bin_size == bin_size_strs[1]) %>% 
  pull(proportion) %>% var
sim_data_child <- sapply(1:length(sim_data), function(k) {
  alpha_k<- sim_data[k]*((sim_data[k]*(1-sim_data[k])/var_prop_child) - 1)
  beta_k <- (1-sim_data[k])*((sim_data[k]*(1-sim_data[k])/var_prop_child) - 1)
  rbeta(n_bins, shape1 = alpha_k, shape2 = beta_k)
}) %>% t() %>% as.data.frame() %>%
  `colnames<-`(paste0("child", 1:4)) %>%
  bind_cols(sim_data %>% data.frame(props_mother=.) %>% mutate(containing_spot=1:n())) %>%
  pivot_longer(cols = starts_with("child"),
               names_to = "child", values_to = "props_child")

sim_data_child_range <- sim_data_child %>%
  group_by(containing_spot) %>%
  mutate(range = max(props_child) - min(props_child))

ggplot(sim_data_child_range %>% ungroup() %>%
         distinct(containing_spot, props_mother, range) %>% 
         arrange(props_mother), 
       aes(x=props_mother, y=range)) +
  geom_abline(slope=1, intercept=0, color="gray70", linetype="dashed") +
  geom_point(size=0.01) +
  theme_classic()

# dir_df <- sapply(1:length(sim_data), function(k) {
#   DirichletReg::rdirichlet(n_bins,
#                            c(sim_data[k], 1-sim_data[k]))[,1]
# }) %>% t() %>% as.data.frame() %>% 
#   `colnames<-`(paste0("child", 1:4)) %>% 
#   bind_cols(sim_data %>% data.frame(props_mother=.) %>% mutate(containing_spot=1:n())) %>% 
#   pivot_longer(cols = starts_with("child"), 
#                names_to = "child", values_to = "props_child")
# 
# DirichletReg::rdirichlet(1000, c(1.15, 1.85)) %>% 
#   data.frame() %>%
#   pivot_longer(cols = everything(), 
#                names_to = "val") %>% 
#   ggplot(aes(x=value, group=val, fill=val)) +
#   geom_histogram(bins=50, position="dodge")
# 
# dir_dist %>% 
#   as.data.frame() %>% 
#   `colnames<-`(paste0("child", 1:4)) %>%
#   pivot_longer(cols = starts_with("child"), 
#                names_to = "child", values_to = "props_child") %>%
#   ggplot(aes(x=child, y=props_child)) +
#   geom_boxplot()



dir_dist <- DirichletReg::rdirichlet(length(sim_data),
                                     (sim_data*n_bins) %>% rep(each=n_bins) %>% 
                                       matrix(ncol=n_bins, byrow=TRUE))

runif(n_bins) < sim_data[1]

dir_df <- dir_dist %>%
  as.data.frame() %>%
  `colnames<-`(paste0("child", 1:4)) %>%
  bind_cols(sim_data %>% data.frame(props_mother=.) %>% mutate(containing_spot=1:n())) %>%
  pivot_longer(cols = starts_with("child"),
               names_to = "child", values_to = "props_child")

dir_df_range <- dir_df %>% 
  group_by(containing_spot) %>% 
  mutate(range = max(props_child) - min(props_child))

ggplot(dir_df_range %>% ungroup() %>% 
         distinct(containing_spot, props_mother, range) %>% 
         arrange(props_mother), 
       aes(x=props_mother, y=range)) +
  geom_abline(slope=1, intercept=0, color="gray70", linetype="dashed") +
  geom_point(size=0.01) +
  theme_classic()


sim_data %>% head
dir_dist %>% head()

dir_df <- sim_data %>% data.frame(proportion=.) %>%
  mutate(dist = "simulated") %>% 
  rbind(deconv_KCs %>% ungroup() %>% filter(bin_size == bin_size_strs[2]) %>% 
          select(proportion) %>% mutate(dist = "deconvolution"))
ggplot(dir_df, aes(x=dist, y= proportion)) +
  geom_boxplot()

mean_prop <- deconv_KCs %>% filter(bin_size == bin_size_strs[1]) %>% 
  pull(proportion) %>% mean
var_prop <- deconv_KCs %>% filter(bin_size == bin_size_strs[1]) %>% 
  pull(proportion) %>% var

alpha <-mean_prop*((mean_prop*(1-mean_prop)/var_prop) - 1)
beta <- (1-mean_prop)*((mean_prop*(1-mean_prop)/var_prop) - 1)

sim_data_child <- rbeta(nrow(mother_bins),
    shape1 = alpha, shape2 = beta)

dir_child_df <- sim_data_child %>% data.frame(proportion=.) %>%
  mutate(dist = "simulated") %>% 
  rbind(mother_bins %>% ungroup() %>% select(prop_child) %>% 
          rename(proportion = prop_child) %>% mutate(dist = "child bins"))

ggplot(dir_child_df, aes(x=dist, y= proportion)) +
  geom_boxplot()

# Randomly link four child bins to one mother bin
sim_data %>% data.frame(props_mother = .) %>% 
  mutate(spot = 1:n())

logit <- function(x){log(x/(1-x))}
sim_data_transformed <- rep(logit(sim_data), each=n_bins)

noise_df <- lapply(seq(0, 1 , by=0.1), function(noise_sd) {
  noise <- rnorm(length(sim_data_transformed), mean=0, sd=noise_sd)
  sim_data_noise <- (1/(1+exp(-(sim_data_transformed + noise)))) %>% 
    data.frame(props_child = .) %>% 
    mutate(containing_spot = rep(1:length(sim_data),
                                 each=n_bins),
           noise = noise_sd) %>% 
    left_join(sim_data %>% data.frame(props_mother = .) %>% 
                mutate(spot = 1:n()),
              by = c("containing_spot" = "spot"))
}) %>% bind_rows()

sim_joined_range <- noise_df %>% 
  group_by(containing_spot, noise) %>% 
  mutate(range = max(props_child) - min(props_child))

ggplot(sim_joined_range %>% ungroup() %>% 
         distinct(containing_spot, props_mother, range, noise) %>% 
         arrange(props_mother), 
       aes(x=props_mother, y=range)) +
  geom_abline(slope=1, intercept=0, color="gray70", linetype="dashed") +
  geom_point(size=0.01) +
  facet_wrap(~noise) +
  theme_classic()

sim_joined <- sim_data_child %>% data.frame(props_child = .) %>% 
  mutate(containing_spot = sample(rep(1:length(sim_data),
                                      each=n_bins))) %>% 
  left_join(sim_data %>% data.frame(props_mother = .) %>% 
              mutate(spot = 1:n()),
            by = c("containing_spot" = "spot"))
# 
# sim_joined <- sim_data_noise %>% 
#   left_join(sim_data %>% data.frame(props_mother = .) %>% 
#               mutate(spot = 1:n()),
#             by = c("containing_spot" = "spot"))

sim_joined_range <- sim_joined %>% 
  group_by(containing_spot) %>% 
  mutate(range = max(props_child) - min(props_child))

ggplot(sim_joined_range %>% ungroup() %>% 
         distinct(containing_spot, props_mother, range) %>% 
         arrange(props_mother), 
       aes(x=props_mother, y=range)) +
  geom_abline(slope=1, intercept=0, color="gray70", linetype="dashed") +
  geom_point(size=0.5)

mother_bins_var <- mother_bins %>% 
  group_by(containing_barcode) %>% 
  # Calculate variance between the four prop_child values and mean = prop_mother
  summarise(var = sum((prop_child - prop_mother)^2) / 3)

top10_variance <- mother_bins_var %>% 
  slice_max(var, n=10)

bottom10_variance <- mother_bins_var %>% 
  slice_min(var, n=10)

mother_bins_top <- mother_bins %>% 
  filter(containing_barcode %in% top10_variance$containing_barcode |
           containing_barcode %in% bottom10_variance$containing_barcode) %>% 
  mutate(containing_barcode = factor(containing_barcode,
                                     levels = c(top10_variance$containing_barcode, 
                                                bottom10_variance$containing_barcode)))

ggplot(mother_bins_top) +
  geom_point(aes(x=as.numeric(containing_barcode), y=prop_child),
             alpha=0.5) +
  geom_linerange(aes(xmin=as.numeric(containing_barcode)-0.5, 
                     xmax=as.numeric(containing_barcode)+0.5,
                     y=prop_mother), color="red")

mother_bins_sample <- mother_bins %>% 
  group_by(containing_barcode) %>% 
  nest() %>%
  ungroup() %>% 
  slice_sample(n=50) %>%
  unnest(cols=c(data)) %>% 
  mutate(containing_barcode = factor(containing_barcode))

ggplot(mother_bins_sample) +
  geom_point(aes(x=as.numeric(containing_barcode), y=prop_child),
             alpha=0.5) +
  geom_linerange(aes(xmin=as.numeric(containing_barcode)-0.5, 
                     xmax=as.numeric(containing_barcode)+0.5,
                     y=prop_mother), color="red")
  
mother_bins_median_var <- mother_bins_var %>% 
  # Get rows where the variance is between 0.05 of median(mother_bins_var$var)
  filter(var > median(var) - 0.0001 & var < median(var) + 0.0001) %>% 
  # Sample 50
  group_by(containing_barcode) %>% 
  nest() %>%
  ungroup() %>% 
  slice_sample(n=50) %>%
  unnest(cols=c(data))

ggplot(mother_bins %>% 
         filter(containing_barcode %in% mother_bins_median_var$containing_barcode) %>% 
         mutate(containing_barcode = factor(containing_barcode))) +
  geom_point(aes(x=as.numeric(containing_barcode), y=prop_child),
             alpha=0.5) +
  geom_linerange(aes(xmin=as.numeric(containing_barcode)-0.5, 
                     xmax=as.numeric(containing_barcode)+0.5,
                     y=prop_mother), color="red") +
  # Plot mean of prop_child
  stat_summary(aes(x=as.numeric(containing_barcode), y=prop_child), 
               fun = mean, geom = "point", color = "blue", size = 2)

# So the mean of all four points isn't equal to the prediction in 16um


tmp <- deconv_props_subset %>% 
  filter(spot %in% mother_bins_median_var$containing_barcode)

ggplot(mother_bins_var, aes(x=var, y = 1)) +
  geom_boxplot()

# histogram
ggplot(mother_bins_var, aes(x=var)) +
  geom_histogram(bins = 100) +
  xlab("Variance of Kupffer cell proportions in mother bins") +
  ylab("Count") +
  ggtitle("Distribution of variance of Kupffer cell proportions in mother bins")

summary(mother_bins_var$var)

