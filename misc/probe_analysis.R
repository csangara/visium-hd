library(tidyverse)

probe_files <- grep(list.files("data/probes", pattern = "*.csv"), pattern = "offtarget",
     invert=TRUE, value=TRUE)

all_probes <- lapply(probe_files, function(file_name){
  offtarget_file <- read.csv(file.path("data/probes", gsub(".csv", ".offtarget.csv", file_name)) )
  probe_file <- read.csv(file.path("data/probes", file_name), skip=5)
  species <- ifelse(grepl("Human", file_name), "Human", "Mouse")
  version <- ifelse(grepl("v2.0", file_name), "2.0", "2.1.0")
  
  probe_file %>% 
    filter(!probe_id %in% offtarget_file$probe_id, included) %>% 
    group_by(gene_id) %>% 
    summarise(num_probes = n()) %>% 
    group_by(num_probes) %>% 
    summarise(num_genes = n()) %>% 
    mutate(percent = (num_genes / sum(num_genes)) * 100) %>% 
    mutate(species = species,
           version = version)
  
}) %>% bind_rows()

all_probes %>% filter(species=="Human") %>% 
  select(num_probes, percent, version) %>% 
  pivot_wider(names_from=version, values_from=percent)

all_probes %>% filter(species=="Mouse") %>% 
  select(num_probes, percent, version) %>% 
  pivot_wider(names_from=version, values_from=percent)

ggplot(all_probes, aes(x=num_probes, y=num_genes, fill=version)) +
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~species) +
  theme_minimal()

all_probes %>% 
  group_by(species, version) %>% 
  summarise(total_genes = sum(num_genes))


# Testing
probes <- read.csv("data/probes/Visium_Human_Transcriptome_Probe_Set_v2.1.0_GRCh38-2024-A.csv", skip=5)
probe_off <- read.csv("data/probes/Visium_Human_Transcriptome_Probe_Set_v2.1.0_GRCh38-2024-A.offtarget.csv")

probes %>% 
  filter(!probe_id %in% probe_off$probe_id) %>% 
  group_by(gene_id) %>% 
  summarise(num_probes = n()) %>% 
  group_by(num_probes) %>% 
  summarise(num_genes = n()) %>% 
  mutate(percent = (num_genes / sum(num_genes)) * 100)
  # ggplot(aes(x=num_probes, y=num_genes)) +
  # geom_bar(stat="identity")

probes$gene_name %>% unique %>% length
probes$gene_id %>% unique %>% length

probes %>% 
  filter(!probe_id %in% probe_off$probe_id,
         included) %>% 
  distinct(gene_id) %>% nrow

genes_more_cov <- probes %>% 
  group_by(gene_id) %>% 
  filter(n() > 3) %>% 
  pull(gene_id)

tmp <- probe_meta %>% 
  filter(gene_id %in% genes_more_cov)

probe_meta %>% 
  tally(included)

probe_meta %>% 
  group_by(transcript_id_set) %>% 
  summarise(num_probes = n()) %>% 
  group_by(num_probes) %>% 
  summarise(num_genes = n()) %>% 
  mutate(percent = num_genes / sum(num_genes) * 100)

probe_meta %>% distinct(gene_name, gene_total_coverage_rounds) %>% 
  group_by(gene_total_coverage_rounds) %>% 
  summarise(num_genes = n()) %>% 
  mutate(percent = num_genes / sum(num_genes) * 100)

probe_genes <- probe_meta$gene_name %>% unique
probe_genes_embl <- probe_meta$gene_id %>% unique

liver_genes <- rownames(liver) 

setdiff(liver_genes, probe_genes)
length(setdiff(liver_genes, probe_genes))
diff_genes_sym <- setdiff(probe_genes, liver_genes)
length(setdiff(probe_genes, liver_genes))

brain <- readRDS("data/Visium_HD_MouseBrain/Visium_HD_MouseBrain_008um.rds")
brain_genes <- rownames(brain)

setdiff(brain_genes, probe_genes_embl)
length(setdiff(brain_genes, probe_genes_embl))
diff_genes_en <- setdiff(probe_genes_embl, brain_genes)
length(setdiff(probe_genes_embl, brain_genes))

probe_meta %>% 
  filter(gene_id %in% diff_genes_en) %>% 
  pull(gene_name) %>% unique %>% 
  setdiff(diff_genes_sym)

liver_caw <- readRDS("data/Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_016um.rds")
liver_caw
setdiff(rownames(liver_caw), rownames(liver))

excluded_probes <- probe_meta %>% 
  group_by(gene_name, gene_description) %>%
  summarise(num_probes = n()) %>% 
  filter(!(gene_name %in% rownames(liver)))
