# Supplementary Figure 4.2
source("visium_hd_liver_combined/0_utils.R")

liver_atlas <- readRDS("~/spotless-benchmark/data/rds/liver_mouseStSt_noEC.rds")

digests <- c("exVivo", "inVivo", "nuclei")
proper_digest_names <- c("scRNA-seq\n(ex vivo digestion)", "scRNA-seq\n(in vivo digestion)", "snRNA-seq") %>%
  setNames(digests)

liver_atlas@meta.data %>%
  # Remove &nbsp; and dots from cell type names
  mutate(annot_cd45 = gsub("[^[:alnum:]]", "", annot_cd45)) %>%
  mutate(annot_cd45 = factor(annot_cd45, levels = rev(celltype_order))) %>% 
  ggplot(aes(y=annot_cd45, x=nCount_RNA, fill=annot_cd45)) +
  geom_boxplot(
    linewidth=0.15,
    outlier.size = 0.4,
    outlier.shape = 16,
    outlier.stroke = 0,
    show.legend = FALSE
  ) +
  stat_summary(fun=mean, geom="point", shape=23, size=0.5, stroke=0.25,
               color="black", fill="white") +
  facet_wrap(~digest,
             scales="free_x",
             labeller = as_labeller(proper_digest_names)
             ) +
  labs(x="Total UMIs per cell") +
  scale_fill_manual(values=color_palette,
                    name="Cell type") +
  scale_y_discrete(labels = celltype_labeller) +
  scale_x_continuous(labels = scales::label_number(scale_cut=scales::cut_short_scale())) +
  theme_minimal(base_size=7) +
  theme_barplot_facet +
  theme(strip.text = element_text(size=7, face="bold"),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        axis.title.x = element_text(size=6,margin=margin(t=5)))

ggsave("visium_hd_liver_combined/plots/boxplot_liver_atlas_counts_per_celltype.pdf",
       width=8, height=5, dpi=300)
