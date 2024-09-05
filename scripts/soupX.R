setwd("/home/chananchidas/topsecret/scripts/")
library(SoupX)

data.dir <- "../spaceranger/binned_outputs/square_008um/"
seuratObj <- load10X(data.dir)

filtered_obj <- Seurat::Read10X(file.path(data.dir, "filtered_feature_bc_matrix"))
raw_obj <- Seurat::Read10X(file.path(data.dir, "raw_feature_bc_matrix"))
raw_obj <- raw_obj[rownames(filtered_obj),]

# Load cluster information by 10x
cluster_info <- read.csv(paste0(data.dir, "analysis_csv/clustering/gene_expression_graphclust/clusters.csv"))
filtered_obj <- filtered_obj[, Cells(filtered_obj) %in% cluster_info$Barcode]

soup_obj <- SoupChannel(raw_obj, filtered_obj)
soup_obj <- setClusters(soup_obj, setNames(as.numeric(cluster_info$Cluster), cluster_info$Barcode))

soup_obj <- autoEstCont(soup_obj)

# To adjust counts
# out <- adjustCounts(soup_obj, roundToInt = TRUE)

# EXTRA: check which genes are not in the filtered object
genes_raw <- rownames(raw_obj) %>% .[!. %in% rownames(filtered_obj)]
dist_df <- raw_obj[genes_raw,] %>% rowSums() %>% data.frame("x"=.)

summary(dist_df)
ggplot(dist_df %>% filter(x < 500), aes(x=x)) + geom_histogram()

