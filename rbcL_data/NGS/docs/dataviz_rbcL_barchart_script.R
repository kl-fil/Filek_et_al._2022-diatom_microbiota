# Dataviz script for Filek et al. 2022 manuscript
# Processing and visualizing data related to rbcL NGS data of source samples
#--------------------load packages--------------------

library(tidyverse) #v1.3.1
library(RColorBrewer)

# set global theme for ggplot
theme_set(theme_light(base_size = 11, base_family = "Lato")) #<< font changes

# define a ggplot theme
theme_barcharts <- theme(
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.9, size = 7, color = "grey30"),
  strip.text.x = element_text(size = 7, color = "gray30", face = "bold"),
  strip.background = element_rect(color = "gray20"),
  panel.grid = element_line(colour = "transparent"),
  panel.grid.major.x = element_line(color = "gray30", linetype = 3),
  panel.border = element_rect(color = "gray30"),
  axis.ticks.x = element_line(colour = "grey30"), # << visibility
  axis.ticks.y = element_line(colour = "grey30"),
  axis.title.x = element_text(size = 8, color = "black"),
  #axis.title.x = element_blank(), #<< visibility
  axis.title.y = element_text(size = 8, color = "black"),
  axis.text.y = element_text(size = 7, color = "grey30"),
  legend.title = element_text(size = 8, color = "black"),
  legend.text = element_text(size = 7, color = "grey30"),
  legend.spacing = unit(0.1, 'cm'),
  legend.key.size = unit(0.3, 'cm'),
  plot.tag = element_text(size = 6)
)

# import taxonomy tables (here I use only lvl 8)

# diat_rbcl_1 <- read.csv("merged_samples/taxonomy/level-1.csv", row.names = 1)
# diat_rbcl_2 <- read.csv("merged_samples/taxonomy/level-2.csv", row.names = 1)
# diat_rbcl_3 <- read.csv("merged_samples/taxonomy/level-3.csv", row.names = 1)
# diat_rbcl_4 <- read.csv("merged_samples/taxonomy/level-4.csv", row.names = 1)
# diat_rbcl_5 <- read.csv("merged_samples/taxonomy/level-5.csv", row.names = 1)
# diat_rbcl_6 <- read.csv("merged_samples/taxonomy/level-6.csv", row.names = 1)
# diat_rbcl_7 <- read.csv("merged_samples/taxonomy/level-7.csv", row.names = 1)
diat_rbcl_8 <- read.csv("merged_samples/taxonomy/level-8.csv", row.names = 1)
# diat_rbcl_9 <- read.csv("merged_samples/taxonomy/level-9.csv", row.names = 1)

metadata_rbcl <- read_tsv("metadata/20210512MonocultbiomeSeqAndSource.tsv") %>%
  rename(SampleID = "#SampleID")

colnames(diat_rbcl_1) #check columns for removal
remove_cols_vect <- colnames(diat_rbcl_1)[2:27]

diat_rbcl_8 %>% #<<
  select(-remove_cols_vect) %>%
  "/"(rowSums(as.data.frame(.))) %>% #calculate proportions
  select(-which(names(.) %in% names(which(apply(., 2, max) < 0.02)))) %>% #select taxa above a percentage / remove taxa below percentage
  select(contains(c("Bacillariophyceae", "Mediophyceae", "Fragilariophyceae", "Coscinodiscophyceae", "Synchromophyceae"))) %>% #select specific taxa #<<
  mutate(Other = 1 - rowSums(.)) %>% #add row with other as 1- sum of abundances
  rename_with(~ gsub("Eukaryota.Chromista.Chromobiota.Bacillariophyta.", "", .x)) %>%
  rename_with(~ gsub("Eukaryota.Sar.Stramenopiles.Ochrophyta.Synchromophyceae.Synchroma.or.Synchroma.fa.Synchroma.ge", "Ochrophyta; Synchromophyceae; Synchroma", .x)) %>%
  rename_with(~ gsub(".__.__.__", " (unclassified)", .x)) %>%
  rename_with(~ gsub(".__", "", .x)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains("TB"), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_rbcl) %>% #<<
  ggplot(aes(x = SampleID, y = abundance)) + 
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", size = 0.35) +
  scale_fill_manual(name = "Diatom taxa (> 2%)", values = c("#A6CEE3", "#1F78B4",
                                                            "#B2DF8A", "#33A02C",
                                                            "#FB9A99", "#E31A1C",
                                                            "#FDBF6F", "#FF7F00",
                                                            "#CAB2D6", "#6A3D9A",
                                                            "#FFFF99", "#B15928",
                                                            "#A6CEE3", "grey70")) +
  scale_y_continuous(labels = function(x) x*100, expand = c(0, 0.01), limits = c(0, 1)) + #set scale size #<<
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)"
  ) +
  theme_barcharts

ggsave("r_output/level-8-taxa-bars.pdf",
       #plot = plotname,
       device = cairo_pdf, #cairo_pdf, png or tiff
       #dpi = 350, #uncomment for png, tiff/ comment for pdf
       width = 180, 
       height = 75, #240 long 150 short
       units = c("mm"))

