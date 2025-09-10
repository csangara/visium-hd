library(Seurat)
library(tidyverse)

liver <- readRDS("../spotless-benchmark/data/rds/liver_mouseStSt_guilliams2022.rds")

# Count plot

ggplot(liver@meta.data, aes(x=nCount_RNA, y=annot_cd45, fill=annot_cd45)) +
  geom_boxplot() +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(~digest) +
  xlab("Counts per spot") +
  ylab("Cell type") +
  ggtitle("Counts per spot by cell type") +
  stat_summary(fun=mean, geom="point", shape=23, size=2, color="black", fill="white")
  