# Dataviz script for Filek et al. 2022 manuscript
## Processing and visualizing data related to initial sequencing results for
## monocultures and source samples combined
##--------------------load packages--------------------

library(tidyverse) #v1.3.1
library(ggplot2) #v3.3.5
library(qiime2R) #qiime to R import qza v0.99.6
library(patchwork) #collect multiple ggplots together v1.1.1
library(ggrepel) #repel text/labels v0.9.1
library(ggforce) #v0.3.3
library(ggtext) #v0.1.1
library(vegan) #v2.5.7
library(ggdendro) #0.1.22


#set global theme for ggplot
theme_set(theme_light(base_size = 11, base_family = "Arial")) #<< font!

# "#<<" comment indicates the line of code that needs to be changed when using 
# a different dataset for the same code

##---- make directory for visualizations ----
dir.create("r_output/visualizations/mono_source_viz", recursive = TRUE)

##---- Load metadata and data-----
# load metadata in tsv format
metadata_ms <- readr::read_tsv("20210512MonocultbiomeSeqAndSource.tsv") %>%
  rename(SampleID = "#SampleID")

# load taxonomies
taxonomy_ms <- read_qza("mono-source/taxonomy/taxonomy.qza")

##---- Load data from qiime2 qza to R for PCoA/PCA plots-----
# load qiime pca and pcoa results.qza (in_data_16 folder)
# 16S data loading (from core-metrics-results-merged_16S-0-with-phyla-no-mitochondria-no-chloroplast-filtered-phylogeny )
# at SeqTry3 folder of 16S analyses
bray_pcoa_ms <- read_qza("mono-source/diversity/core-metrics-results-D32260/bray_curtis_pcoa_results.qza")
unifrac_pcoa_ms <- read_qza("mono-source/diversity/core-metrics-results-D32260/unweighted_unifrac_pcoa_results.qza")
w_unifrac_pcoa_ms <- read_qza("mono-source/diversity/core-metrics-results-D32260/weighted_unifrac_pcoa_results.qza")
jaccard_pcoa_ms <- read_qza("mono-source/diversity/core-metrics-results-D32260/jaccard_pcoa_results.qza")
rAitch_pca_ms <- read_qza("mono-source/diversity/deicode/rAitchison-ordination.qza")
g_unifrac_pcoa_ms <- read_qza("mono-source/diversity/gen-unifrac-D32260/generalized_unifrac_pcoa_results.qza")

##---- PCoA plots ----

# define theme for pcoa/pca plots
theme_pcoa1 <- theme(
  panel.border = element_rect(fill = "transparent", color = "grey30"),
  plot.subtitle = element_text(size = 9),
  axis.title.x = element_text(size = 9, color = "black"),
  axis.title.y = element_text(size = 9, color = "black"),
  axis.text.x = element_text(size = 7, color = "grey30"),
  axis.text.y = element_text(size = 7, color = "grey30"),
  legend.title = element_text(size = 9, color = "black"),
  legend.text = element_markdown(size = 8, color = "grey30"),
  legend.spacing.y = unit(0.05, 'cm'),
  legend.key.size = unit(0.3, 'cm'),
  legend.position = "bottom",
  legend.box = "horizontal",
  panel.grid.minor.x = element_line(size = 0, color = "transparent"),
  panel.grid.minor.y = element_line(size = 0, color = "transparent"),
  panel.grid.major.x = element_line(size = 0, color = "transparent"),
  panel.grid.major.y = element_line(size = 0, color = "transparent"),
  axis.line.x.bottom = element_line(size = 0.3, color = "grey30"),
  axis.line.y.left = element_line(size = 0.3, color = "grey30"),
  axis.ticks = element_line(colour = "grey30")
)

##---- bray1_ms----
bray1_plot_ms <-
  bray_pcoa_ms$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_ms) %>%  #<<
  mutate(Genus = fct_relevel(Genus,
                             "Achnanthes", "Amphora", "Poulinea",
                             "Diploneis", "Fallacia", "Nitzschia",
                             "Entomoneis", "Psammodictyon", "Caretta")) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = Genus, shape = TurtleID),
    alpha=0.8,
    size = 3
  ) +
  geom_point(
    aes(shape = TurtleID),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Turtle ID:") +
  scale_fill_manual(name = "Genus:", 
                    values = c("#1B9E77", "#D95F02", "#7570B3", 
                               "#E7298A", "#66A61E", "#E6AB02", 
                               "#A6761D", "#666666", "#FFFFFF"),
                    labels = c(Achnanthes = "*Achnanthes*",
                               Amphora = "*Amphora*",
                               Poulinea = "*Poulinea*",
                               Diploneis = "*Diploneis*",
                               Fallacia = "*Fallacia*",
                               Nitzschia = "*Nitzschia*",
                               Entomoneis = "*Entomoneis*",
                               Psammodictyon = "*Psammodictyon*",
                               Caretta = "*Caretta* (host)")) +
  scale_y_continuous(breaks = c(-0.4, -0.2, 0, 0.2)) +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2, 0.4)) +
  guides(fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2),
         shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1)) +
  xlab(paste("PC1 (", round(100*bray_pcoa_ms$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*bray_pcoa_ms$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Bray-Curtis")

##---- jacc1_ms----
jacc1_plot_ms <-
  jaccard_pcoa_ms$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_ms) %>%  #<<
  mutate(Genus = fct_relevel(Genus,
                             "Achnanthes", "Amphora", "Poulinea",
                             "Diploneis", "Fallacia", "Nitzschia",
                             "Entomoneis", "Psammodictyon", "Caretta")) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = Genus, shape = TurtleID),
    alpha=0.8,
    size = 3
  ) +
  geom_point(
    aes(shape = TurtleID),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Turtle ID:") +
  scale_fill_manual(name = "Genus:", 
                    values = c("#1B9E77", "#D95F02", "#7570B3", 
                               "#E7298A", "#66A61E", "#E6AB02", 
                               "#A6761D", "#666666", "#FFFFFF"),
                    labels = c(Achnanthes = "*Achnanthes*",
                               Amphora = "*Amphora*",
                               Poulinea = "*Poulinea*",
                               Diploneis = "*Diploneis*",
                               Fallacia = "*Fallacia*",
                               Nitzschia = "*Nitzschia*",
                               Entomoneis = "*Entomoneis*",
                               Psammodictyon = "*Psammodictyon*",
                               Caretta = "*Caretta* (host)")) +
  scale_y_continuous(breaks = c(-0.2, 0, 0.2)) +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2, 0.4)) +
  guides(fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2),
         shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1)) +
  xlab(paste("PC1 (", round(100*jaccard_pcoa_ms$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*jaccard_pcoa_ms$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Jaccard")

##---- wuni1_ms----
wuni1_plot_ms <-
  w_unifrac_pcoa_ms$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_ms) %>%  #<<
  mutate(Genus = fct_relevel(Genus,
                             "Achnanthes", "Amphora", "Poulinea",
                             "Diploneis", "Fallacia", "Nitzschia",
                             "Entomoneis", "Psammodictyon", "Caretta")) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = Genus, shape = TurtleID),
    alpha=0.8,
    size = 3
  ) +
  geom_point(
    aes(shape = TurtleID),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Turtle ID:") +
  scale_fill_manual(name = "Genus:", 
                    values = c("#1B9E77", "#D95F02", "#7570B3", 
                               "#E7298A", "#66A61E", "#E6AB02", 
                               "#A6761D", "#666666", "#FFFFFF"),
                    labels = c(Achnanthes = "*Achnanthes*",
                               Amphora = "*Amphora*",
                               Poulinea = "*Poulinea*",
                               Diploneis = "*Diploneis*",
                               Fallacia = "*Fallacia*",
                               Nitzschia = "*Nitzschia*",
                               Entomoneis = "*Entomoneis*",
                               Psammodictyon = "*Psammodictyon*",
                               Caretta = "*Caretta* (host)")) +
  scale_y_continuous(breaks = c(-0.2, -0.1, 0, 0.1)) +
  scale_x_continuous(breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3)) +
  guides(fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2),
         shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1)) +
  xlab(paste("PC1 (", round(100*w_unifrac_pcoa_ms$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*w_unifrac_pcoa_ms$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Weighted Unifrac")

##---- uni1_ms----
uni1_plot_ms <-
  unifrac_pcoa_ms$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_ms) %>%  #<<
  mutate(Genus = fct_relevel(Genus,
                             "Achnanthes", "Amphora", "Poulinea",
                             "Diploneis", "Fallacia", "Nitzschia",
                             "Entomoneis", "Psammodictyon", "Caretta")) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = Genus, shape = TurtleID),
    alpha=0.8,
    size = 3
  ) +
  geom_point(
    aes(shape = TurtleID),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Turtle ID:") +
  scale_fill_manual(name = "Genus:", 
                    values = c("#1B9E77", "#D95F02", "#7570B3", 
                               "#E7298A", "#66A61E", "#E6AB02", 
                               "#A6761D", "#666666", "#FFFFFF"),
                    labels = c(Achnanthes = "*Achnanthes*",
                               Amphora = "*Amphora*",
                               Poulinea = "*Poulinea*",
                               Diploneis = "*Diploneis*",
                               Fallacia = "*Fallacia*",
                               Nitzschia = "*Nitzschia*",
                               Entomoneis = "*Entomoneis*",
                               Psammodictyon = "*Psammodictyon*",
                               Caretta = "*Caretta* (host)"))+
  guides(fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2),
         shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1)) +
  xlab(paste("PC1 (", round(100*unifrac_pcoa_ms$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*unifrac_pcoa_ms$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Unweighted Unifrac")

##---- guni1_ms----
guni1_plot_ms <-
  g_unifrac_pcoa_ms$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_ms) %>%  #<<
  mutate(Genus = fct_relevel(Genus,
                             "Achnanthes", "Amphora", "Poulinea",
                             "Diploneis", "Fallacia", "Nitzschia",
                             "Entomoneis", "Psammodictyon", "Caretta")) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = Genus, shape = TurtleID),
    alpha=0.8,
    size = 3
  ) +
  geom_point(
    aes(shape = TurtleID),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Turtle ID:") +
  scale_fill_manual(name = "Genus:", 
                    values = c("#1B9E77", "#D95F02", "#7570B3", 
                               "#E7298A", "#66A61E", "#E6AB02", 
                               "#A6761D", "#666666", "#FFFFFF"),
                    labels = c(Achnanthes = "*Achnanthes*",
                               Amphora = "*Amphora*",
                               Poulinea = "*Poulinea*",
                               Diploneis = "*Diploneis*",
                               Fallacia = "*Fallacia*",
                               Nitzschia = "*Nitzschia*",
                               Entomoneis = "*Entomoneis*",
                               Psammodictyon = "*Psammodictyon*",
                               Caretta = "*Caretta* (host)"))+
  guides(fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2),
         shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1)) +
  xlab(paste("PC1 (", round(100*g_unifrac_pcoa_ms$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*g_unifrac_pcoa_ms$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Generalized Unifrac")

##---- bray2_ms----
bray2_plot_ms <-
  bray_pcoa_ms$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_ms) %>%  #<<
  mutate(Genus = fct_relevel(Genus,
                             "Achnanthes", "Amphora", "Poulinea",
                             "Diploneis", "Fallacia", "Nitzschia",
                             "Entomoneis", "Psammodictyon", "Caretta")) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = SourceIDfig, shape = TypeExp),
    alpha=0.8,
    size = 3
  ) +
  geom_point(
    aes (shape = TypeExp),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  scale_fill_manual(name = "Source sample ID:", values = c("#1F78B4", 
                                                           "#33A02C", "#B2DF8A", 
                                                           "#FF7F00", 
                                                           "#6A3D9A", "#CAB2D6")) +
  scale_shape_manual(values = c(21, 24), name = "Sample type:") +
  guides(fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2),
         shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1)) +
  scale_y_continuous(breaks = c(-0.4, -0.2, 0, 0.2)) +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2, 0.4)) +
  xlab(paste("PC1 (", round(100*bray_pcoa_ms$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*bray_pcoa_ms$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #set theme
  labs(subtitle = "Bray-Curtis")

##---- jacc2_ms----
jacc2_plot_ms <-
  jaccard_pcoa_ms$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_ms) %>%  #<<
  mutate(Genus = fct_relevel(Genus,
                             "Achnanthes", "Amphora", "Poulinea",
                             "Diploneis", "Fallacia", "Nitzschia",
                             "Entomoneis", "Psammodictyon", "Caretta")) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = SourceIDfig, shape = TypeExp),
    alpha=0.8,
    size = 3
  ) +
  geom_point(
    aes (shape = TypeExp),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  scale_fill_manual(name = "Source sample ID:", values = c("#1F78B4", 
                                                           "#33A02C", "#B2DF8A", 
                                                           "#FF7F00", 
                                                           "#6A3D9A", "#CAB2D6")) +
  scale_shape_manual(values = c(21, 24), name = "Sample type:") +
  guides(fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2),
         shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1)) +
  xlab(paste("PC1 (", round(100*jaccard_pcoa_ms$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*jaccard_pcoa_ms$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #set theme
  labs(subtitle = "Jaccard")

##---- wuni2_ms----
wuni2_plot_ms <-
  w_unifrac_pcoa_ms$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_ms) %>%  #<<
  mutate(Genus = fct_relevel(Genus,
                             "Achnanthes", "Amphora", "Poulinea",
                             "Diploneis", "Fallacia", "Nitzschia",
                             "Entomoneis", "Psammodictyon", "Caretta")) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = SourceIDfig, shape = TypeExp),
    alpha=0.8,
    size = 3
  ) +
  geom_point(
    aes (shape = TypeExp),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  scale_fill_manual(name = "Source sample ID:", values = c("#1F78B4", 
                                                           "#33A02C", "#B2DF8A", 
                                                           "#FF7F00", 
                                                           "#6A3D9A", "#CAB2D6")) +
  scale_shape_manual(values = c(21, 24), name = "Sample type:") +
  guides(fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2),
         shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1)) +
  xlab(paste("PC1 (", round(100*w_unifrac_pcoa_ms$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*w_unifrac_pcoa_ms$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #set theme
  labs(subtitle = "Weighted Unifrac")

##---- uni2_ms----
uni2_plot_ms <-
  unifrac_pcoa_ms$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_ms) %>%  #<<
  mutate(Genus = fct_relevel(Genus,
                             "Achnanthes", "Amphora", "Poulinea",
                             "Diploneis", "Fallacia", "Nitzschia",
                             "Entomoneis", "Psammodictyon", "Caretta")) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = SourceIDfig, shape = TypeExp),
    alpha=0.8,
    size = 3
  ) +
  geom_point(
    aes (shape = TypeExp),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  scale_fill_manual(name = "Source sample ID:", values = c("#1F78B4", 
                                                           "#33A02C", "#B2DF8A", 
                                                           "#FF7F00", 
                                                           "#6A3D9A", "#CAB2D6")) +
  scale_shape_manual(values = c(21, 24), name = "Sample type:") +
  guides(fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2),
         shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1)) +
  xlab(paste("PC1 (", round(100*unifrac_pcoa_ms$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*unifrac_pcoa_ms$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #set theme
  labs(subtitle = "Unweighted Unifrac")

##---- guni2_ms----
guni2_plot_ms <-
  g_unifrac_pcoa_ms$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_ms) %>%  #<<
  mutate(Genus = fct_relevel(Genus,
                             "Achnanthes", "Amphora", "Poulinea",
                             "Diploneis", "Fallacia", "Nitzschia",
                             "Entomoneis", "Psammodictyon", "Caretta")) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = SourceIDfig, shape = TypeExp),
    alpha=0.8,
    size = 3
  ) +
  geom_point(
    aes (shape = TypeExp),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  scale_fill_manual(name = "Source sample ID:", values = c("#1F78B4", 
                                                           "#33A02C", "#B2DF8A", 
                                                           "#FF7F00", 
                                                           "#6A3D9A", "#CAB2D6")) +
  scale_shape_manual(values = c(21, 24), name = "Sample type:") +
  guides(fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2),
         shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1)) +
  xlab(paste("PC1 (", round(100*g_unifrac_pcoa_ms$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*g_unifrac_pcoa_ms$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1  + #set theme
  labs(subtitle = "Generalized Unifrac")

##---- PCA plots-----
# rename column in taxonomy so we can join it to our data
taxonomy_ms <- rename(taxonomy_ms$data, FeatureID = Feature.ID)

rPCA1_plot_ms <- 
  ggplot() +
  geom_segment(data = rAitch_pca_ms$data$Species %>% #<<
                 mutate(a = sqrt(PC1^2+PC2^2+PC3^2)) %>% #calculate the distance from the origin
                 slice_max(n = 5, order_by = a) %>% #keep X furthest away points - instead of top_n -- slice_max
                 mutate(PC1=PC1*4.5, PC2=PC2*4.5) %>% #scale arrow linearly, i guess it's ok (look at emperor)
                 left_join(taxonomy_ms), #<<
               aes(x=0, xend=PC1, y=0, yend=PC2),
               size = 0.3,
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray30") +
  geom_point(
    data = rAitch_pca_ms$data$Vectors %>% #<<
      left_join(metadata_ms) %>%
      mutate(Genus = fct_relevel(Genus,
                                 "Achnanthes", "Amphora", "Poulinea",
                                 "Diploneis", "Fallacia", "Nitzschia",
                                 "Entomoneis", "Psammodictyon", "Caretta")), #<<
    aes(
      x = PC1, 
      y = PC2, 
      fill = Genus, 
      shape = TurtleID
    ),
    alpha = 0.8,
    size = 3
  ) +
  geom_point(
    data = rAitch_pca_ms$data$Vectors %>% #<<
      left_join(metadata_ms), #<<
    aes(
      x = PC1, 
      y = PC2, 
      shape = TurtleID
    ),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  geom_text_repel(data = rAitch_pca_ms$data$Species %>% #<<
                    mutate(a = sqrt(PC1^2+PC2^2+PC3^2)) %>% #calculate the distance from the origin
                    slice_max(n = 5, order_by = a) %>% #keep X furthest away points - instead of top_n -- slice_max
                    mutate(PC1=PC1*4.5, PC2=PC2*4.5) %>% #scale arrow linearly, i guess it's ok (look at emperor)
                    left_join(taxonomy_ms) %>% #<<
                    mutate(Taxa = c("plain('uncultured')~plain('Oligoflexales')",
                                    "plain('')~italic('Flavobacteriaceae')",
                                    "plain('')~italic('Rhodobacteraceae')",
                                    "plain('')~italic('Rhodobacteraceae')",
                                    "plain('')~italic('Rhodobacteraceae')")),
                  aes(x = PC1, y = PC2, 
                      label = Taxa),
                  point.padding = 2,
                  alpha = 0.8,
                  colour = "black",
                  size = 2,
                  parse = TRUE) +
  scale_shape_manual(name = "Turtle ID:", values = c(21, 22, 23, 24)) +
  scale_fill_manual(name = "Genus:", 
                    values = c("#1B9E77", "#D95F02", "#7570B3", 
                               "#E7298A", "#66A61E", "#E6AB02", 
                               "#A6761D", "#666666", "#FFFFFF"),
                    labels = c(Achnanthes = "*Achnanthes*",
                               Amphora = "*Amphora*",
                               Poulinea = "*Poulinea*",
                               Diploneis = "*Diploneis*",
                               Fallacia = "*Fallacia*",
                               Nitzschia = "*Nitzschia*",
                               Entomoneis = "*Entomoneis*",
                               Psammodictyon = "*Psammodictyon*",
                               Caretta = "*Caretta* (host)")) +
  guides(shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1), 
         fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2)) +
  xlab(paste("PC1 (",round(100*rAitch_pca_ms$data$ProportionExplained[1],2),"%)")) + #<<
  ylab(paste("PC2 (",round(100*rAitch_pca_ms$data$ProportionExplained[2],2),"%)")) + #<<
  theme_pcoa1 +
  labs(subtitle = "Robust Aitchison")

rPCA2_plot_ms <-
  rAitch_pca_ms$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_ms) %>%  #<<
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = SourceIDfig, shape = TypeExp),
    alpha=0.8,
    size = 3
  ) +
  geom_point(
    aes (shape = TypeExp),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  scale_fill_manual(name = "Source sample ID:", values = c("#1F78B4", 
                                                           "#33A02C", "#B2DF8A", 
                                                           "#FF7F00", 
                                                           "#6A3D9A", "#CAB2D6")) +
  scale_shape_manual(values = c(21, 24), name = "Sample type:") +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2, 0.4)) +
  guides(fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2),
         shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1)) +
  xlab(paste("PC1 (", round(100*rAitch_pca_ms$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*rAitch_pca_ms$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 +
  labs(subtitle = "Robust Aitchison")



##---- saving pcoa/pca plots -----

diss_genus_plot <- (bray1_plot_ms + wuni1_plot_ms) / (jacc1_plot_ms + uni1_plot_ms) /
  (rPCA1_plot_ms + guni1_plot_ms) +
  plot_annotation(tag_levels = 'A', tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect")

diss_source_plot <- (bray2_plot_ms + wuni2_plot_ms) / (jacc2_plot_ms + uni2_plot_ms) /
  (rPCA2_plot_ms + guni2_plot_ms) +
  plot_annotation(tag_levels = 'A', tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect")

ggsave("r_output/visualizations/mono_source_viz/suppl_dissim_all_genus.pdf",
       plot = diss_genus_plot,
       #dpi = ,
       device = cairo_pdf,
       width = 168,
       height = 200, 
       units = c("mm"))

ggsave("r_output/visualizations/mono_source_viz/suppl_dissim_all_source.pdf",
       plot = diss_source_plot,
       #dpi = ,
       device = cairo_pdf,
       width = 168,
       height = 200, 
       units = c("mm"))

#add margins to the right of the plots (for 1st column in multipanel plots)

bray1_plot_ms_a <- bray1_plot_ms +
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))

jacc1_plot_ms_a <- jacc1_plot_ms +
  #   subtitle = "Colored by source sample"
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))

wuni1_plot_ms_a <- wuni1_plot_ms +
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))

uni1_plot_ms_a <- uni1_plot_ms +
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))

guni1_plot_ms_a <- guni1_plot_ms +
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))

rPCA1_plot_ms_a <- rPCA1_plot_ms +
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))

all_ord_plots2_ms_a <- (bray1_plot_ms_a + bray2_plot_ms + jacc1_plot_ms_a + jacc2_plot_ms) /
  (wuni1_plot_ms_a + wuni2_plot_ms + uni1_plot_ms_a + uni2_plot_ms) /
  (guni1_plot_ms_a + guni2_plot_ms + rPCA1_plot_ms_a + rPCA2_plot_ms) +
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        plot.tag.position = c(0.015,0.995),
        plot.tag = element_text(size = 10))

ggsave("r_output/visualizations/mono_source_viz/all-250-600-test.pdf",
       plot = all_ord_plots2_ms_a,
       device = cairo_pdf, #cairo_pdf, png or tiff
       #dpi = 350, #uncomment for png, tiff/ comment for pdf
       width = 250, 
       height = 600, #240 long 150 short
       units = c("mm"))

##---- Dendrograms - generalized unifrac ----
# import distance matrices (generalized unifrac)
guni_dist_ms <- read_qza("mono-source/diversity/gen-unifrac-D32260/generalized_unifrac_distance_matrix.qza")

# plot generalized unifrac clustering
plot(hclust(guni_dist_ms$data, method = "average"), hang = -1, cex = 1) 

guni_d_ms <- hclust(guni_dist_ms$data, method = "average")

# extract order of samples
guni_d_order_labels_ms <- guni_d_ms$labels[guni_d_ms$order]
guni_d_order_labels_ms_dm <- guni_d_order_labels_ms %>%
  str_replace("KF", "DM")

ggdendrogram(guni_d_ms) #check

# turn to dendrogram for ggplot
d_guni_d_ms <- as.dendrogram(guni_d_ms) %>%
  dendro_data(., type = "rectangle")

# test visualization
ggplot(segment(d_guni_d_ms)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))

# produce dendrogram plot with sample symbols
p_dendro1 <- 
  d_guni_d_ms$segments %>% #<<
  mutate(yend2 = case_when(yend == 0 ~ yend + 0.2,
                           yend != 0 ~ yend)) %>%
  ggplot() + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend2), lineend = "round") +
  scale_x_continuous(expand = c(0.027,0)) + #<< handle positioning here to align with bar charts
  scale_y_continuous(expand = c(0.12,0)) + #<<
  geom_point(data = metadata_ms %>% #<<
               mutate(label = fct_relevel(SampleID, guni_d_order_labels_ms)) %>% #<<
               left_join(d_guni_d_ms$labels, by = "label") %>%
               mutate(y2 = case_when(y == 0 ~ y + 0.2)) %>%
               mutate(Genus = fct_relevel(Genus,
                                          "Achnanthes", "Amphora", "Poulinea",
                                          "Diploneis", "Fallacia", "Nitzschia",
                                          "Entomoneis", "Psammodictyon", "Caretta")),
             aes(x = x, y = y2, fill = Genus, shape = TurtleID),
             alpha = 1,
             size = 3) +
  scale_shape_manual(name = "Turtle ID:", values = c(21, 22, 23, 24)) +
  scale_fill_manual(name = "Genus:", 
                    values = c("#1B9E77", "#D95F02", "#7570B3", 
                               "#E7298A", "#66A61E", "#E6AB02", 
                               "#A6761D", "#666666", "#FFFFFF"),
                    labels = c(Achnanthes = "*Achnanthes*",
                               Amphora = "*Amphora*",
                               Poulinea = "*Poulinea*",
                               Diploneis = "*Diploneis*",
                               Fallacia = "*Fallacia*",
                               Nitzschia = "*Nitzschia*",
                               Entomoneis = "*Entomoneis*",
                               Psammodictyon = "*Psammodictyon*",
                               Caretta = "*Caretta* (host)")) +
  guides(shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1), 
         fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2)) +
  labs(tag = "(A)") +
  theme_void() + 
  theme(legend.text = element_markdown(size = 8, color = "grey30"),
        legend.position = "none", # show or do not show legend
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0.1,0), "cm"),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key.size = unit(0.3, 'cm'),
        plot.tag = element_text(size = 12)) +
  coord_fixed(ratio = 10, expand = TRUE, clip = "on")

# plot with source sample colorings  
p_dendro2 <-
  ggplot() +  
  geom_point(data = metadata_ms %>%
               mutate(label = fct_relevel(SampleID, guni_d_order_labels_ms)) %>%
               left_join(d_guni_d_ms$labels, by = "label") %>%
               mutate(y2 = case_when(y == 0 ~ y + 0.175)),
             aes(x = x, y = y2, fill = SourceIDfig),
             alpha = 1,
             size = 3,
             shape = 21) +
  scale_fill_manual(name = "Source sample ID:", values = c("#1F78B4", 
                                                           "#33A02C", "#B2DF8A",
                                                           "#FF7F00",
                                                           "#6A3D9A", "#CAB2D6")) +
  scale_x_continuous(expand = c(0.027,0)) +
  scale_y_continuous(expand = c(0.3,0)) +
  theme_void() +
  theme(legend.text = element_markdown(size = 8, color = "grey30"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0.1,0), "cm"),
        plot.tag = element_text(size = 6)) +
  coord_fixed()

# combine plots
p_dendro_combination <- p_dendro1/p_dendro2

p_dendro_combination <- p_dendro_combination +
  coord_fixed(ratio = 1, expand = TRUE, clip = "on")





##---- Import taxonomy and counts (extracted from qiime2 barplots) ----

# import csv files taxa as columns, samples as rows
data_2_ms <- read.csv("mono-source/taxonomy/level-2.csv", row.names = 1) #row names as sample id
data_3_ms <- read.csv("mono-source/taxonomy/level-3.csv", row.names = 1) 
data_4_ms <- read.csv("mono-source/taxonomy/level-4.csv", row.names = 1) 
data_5_ms <- read.csv("mono-source/taxonomy/level-5.csv", row.names = 1)
data_6_ms <- read.csv("mono-source/taxonomy/level-6.csv", row.names = 1)

# extract taxa names present only in diatom monocultures

data_2_ms[1:19,] %>% #select first 19 rows (monoculture samples)
  "/"(rowSums(.)) #divide by row sums to get relative abundance

data_2_ms[1:19,] %>% #select first 19 rows (monoculture samples)
  colSums() #sum taxa reads over all samples

# select taxa names present only in monocultures in data_X_ms level
mononames_2 <- names(
  which(data_2_ms[1:19,] %>%
          apply(2,max) > 0)
)

# normalize counts - average per column - proportions
data_2_ms_prop <- data_2_ms/rowSums(data_2_ms)
data_3_ms_prop <- data_3_ms/rowSums(data_3_ms)
data_4_ms_prop <- data_4_ms/rowSums(data_4_ms)
data_5_ms_prop <- data_5_ms/rowSums(data_5_ms)
data_6_ms_prop <- data_6_ms/rowSums(data_6_ms)

# select taxa names present only in monocultures in data_X_ms level
mononames_2 <- names(
  which(data_2_ms[1:19,] %>%
          apply(2,max) > 0)
)

data_2_ms_prop %>%
  select(mononames_2)

# same thing as above without assigning external variables
data_2_ms_prop %>% #<<
  select(names(
    which(data_2_ms[1:19,] %>%
            apply(2,max) > 0)
  ))

##---- Barplots - monocultures and source samples ----

# define general theme to use for plots
theme_barcharts <- theme(
  axis.text.x = element_text(angle = 90, vjust = 0.5, 
                             hjust = 0.9, size = 11, 
                             color = "grey30"), 
  strip.text.x = element_text(size = 11, color = "gray30", face = "bold"),
  strip.background = element_rect(color = "gray20"),
  panel.grid = element_line(colour = "transparent"),
  panel.grid.major.x = element_line(color = "gray30", linetype = 3),
  panel.border = element_rect(color = "gray30"),
  axis.ticks.x = element_line(colour = "grey30"), # << visibility
  axis.ticks.y = element_line(colour = "grey30"),
  axis.title.x = element_text(size = 11, color = "black"),
  #axis.title.x = element_blank(), #<< visibility
  axis.title.y = element_text(size = 11, color = "black"),
  axis.text.y = element_text(size = 9, color = "grey30"),
  legend.title = element_text(size = 9, color = "black"),
  legend.text = element_text(size = 8, color = "grey30"),
  legend.spacing.y = unit(0.05, 'cm'),
  legend.key.size = unit(0.3, 'cm'),
  plot.tag = element_text(size = 12)
)

taxaplot_2 <- 
  data_2_ms %>% #<<
  "/"(rowSums(.)) %>% #calculate proportions
  select(names( #select names from
    which(.[1:19,] %>% #monoculture samples only
            apply(2,max) > 0) #above max value zero (select taxa with counts)
  )) %>%
  select(-which(names(.) %in% names(which(apply(., 2, max) < 0.01)))) %>% #select taxa above a percentage / remove taxa below percentage
  mutate(Other = 1 - rowSums(.)) %>% #add row with other as 1- sum of abundances
  rename_with(~ gsub("d__Bacteria.p__", " ", .x)) %>%
  t() %>% #transpose
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("TB", "KF")), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_ms) %>% #<<
  mutate(SamplePlotName = fct_relevel(SamplePlotName, guni_d_order_labels_ms_dm)) %>%
  mutate(ExpGen = paste(TypeExp, Genus, sep = "_")) %>%
  ggplot(aes(x = SamplePlotName, y = abundance)) + 
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", size = 0.35) +
  scale_fill_brewer(name = "Phyla (> 1%)", palette = "Paired") +
  scale_y_continuous(breaks = c(0, 0.25, .50, .75, 1), 
                     labels = function(x) x*100, expand = c(0, 0), 
                     limits = c(0, 1)) + #set scale size #<<
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)"
  ) +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme_barcharts +
  coord_fixed(ratio = 10, expand = TRUE, clip = "on")

taxaplot_2_nox <- taxaplot_2 +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(tag = "(B)")

taxaplot_alpha_4 <- 
  data_5_ms %>% #<<
  "/"(rowSums(.)) %>% #calculate proportions
  select(names( #select names from
    which(.[1:19,] %>% #monoculture samples only
            apply(2,max) > 0) #above max value zero (select taxa with counts)
  )) %>%
  select(-which(names(.) %in% names(which(apply(., 2, max) < 0.05)))) %>% #select taxa above a percentage / remove taxa below percentage
  select(contains("Alphaproteobacteria")) %>% #select specific taxa #<<
  rename_with(~ gsub("d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria", " ", .x)) %>%
  rename_with(~ gsub(".o__", "", .x)) %>%
  rename_with(~ gsub(".f__", "; ", .x)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("TB", "KF")), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_ms) %>% #<<
  mutate(SamplePlotName = fct_relevel(SamplePlotName, guni_d_order_labels_ms_dm)) %>%
  mutate(ExpGen = paste(TypeExp, Genus, sep = "_")) %>%
  ggplot(aes(x = SamplePlotName, y = abundance)) +  
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", size = 0.35) +
  scale_fill_manual(name = "Alphaproteobacteria (> 5%)", values = c("#A6CEE3", "#1F78B4",                                                                                #                                                                    "#B2DF8A", "#33A02C",
                                                                    "#FB9A99", "#E31A1C", 
                                                                    "#FDBF6F", "#FF7F00",
                                                                    "#CAB2D6", "#6A3D9A",
                                                                    "#FFFF99", "#B15928",
                                                                    "#A6CEE3", "#1F78B4",
                                                                    "#B2DF8A", "#E31A1C", 
                                                                    "#FDBF6F")) +
  scale_y_continuous(labels = function(x) x*100, expand = c(0, 0), 
                     limits = c(0, 0.45)) + #set scale size #<<
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)"
  ) +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme_barcharts +
  coord_fixed(ratio = 22, expand = TRUE, clip = "on")

taxaplot_alpha_4_nox <- taxaplot_alpha_4 +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank()
  ) +
  labs(tag = "(C)")

taxaplot_gamma_4 <- 
  data_5_ms %>% #<<
  "/"(rowSums(.)) %>% #calculate proportions
  select(names( #select names from
    which(.[1:19,] %>% #monoculture samples only
            apply(2,max) > 0) #above max value zero (select taxa with counts)
  )) %>%
  select(-which(names(.) %in% names(which(apply(., 2, max) < 0.05)))) %>% #select taxa above a percentage / remove taxa below percentage
  select(contains("Gammaproteobacteria")) %>% #select specific taxa #<<
  rename_with(~ gsub("d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria", "", .x)) %>%
  rename_with(~ gsub(".__.__", "Unclassified", .x)) %>%
  rename_with(~ gsub(".o__", "", .x)) %>%
  rename_with(~ gsub(".f__", " ", .x)) %>%
  rename_with(~ gsub(" ", "; ", .x)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("TB", "KF")), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_ms) %>% #<<
  mutate(SamplePlotName = fct_relevel(SamplePlotName, guni_d_order_labels_ms_dm)) %>%
  mutate(ExpGen = paste(TypeExp, Genus, sep = "_")) %>%
  ggplot(aes(x = SamplePlotName, y = abundance)) +  
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", size = 0.35) +
  scale_fill_brewer(name = "Gammaproteobacteria (> 5%)", palette = "Paired") +
  scale_y_continuous(breaks = c(0, 0.25, .50, .75), 
                     labels = function(x) x*100, 
                     expand = c(0, 0), limits = c(0, 1)) + #set scale size #<<
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)"
  ) +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme_barcharts +
  coord_fixed(ratio = 10, expand = TRUE, clip = "on")

taxaplot_gamma_4_nox <- taxaplot_gamma_4 +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.title.y = element_blank()
  ) +
  labs(tag = "(D)")

taxaplot_bacter_4 <- 
  data_5_ms %>% #<<
  "/"(rowSums(.)) %>% #calculate proportions
  select(names( #select names from
    which(.[1:19,] %>% #monoculture samples only
            apply(2,max) > 0) #above max value zero (select taxa with counts)
  )) %>%
  select(-which(names(.) %in% names(which(apply(., 2, max) < 0.01)))) %>% #select taxa above a percentage / remove taxa below percentage
  select(contains("Bacteroidota")) %>% #select specific taxa #<<
  rename_with(~ gsub("d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__", "", .x)) %>%
  rename_with(~ gsub("d__Bacteria.p__Bacteroidota.c__Kapabacteria.o__", "", .x)) %>%
  rename_with(~ gsub("d__Bacteria.p__Bacteroidota.c__Rhodothermia.o__", "", .x)) %>%
  rename_with(~ gsub(".f__", "; ", .x)) %>%
  rename_with(~ gsub("_marine_group", "", .x)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("TB", "KF")), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_ms) %>% #<<
  mutate(SamplePlotName = fct_relevel(SamplePlotName, guni_d_order_labels_ms_dm)) %>%
  mutate(ExpGen = paste(TypeExp, Genus, sep = "_")) %>%
  ggplot(aes(x = SamplePlotName, y = abundance)) +  
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", size = 0.35) +
  scale_fill_brewer(name = "Phylum Bacteroidota (> 1%)", palette = "Paired") +
  scale_y_continuous(labels = function(x) x*100, 
                     expand = c(0, 0), 
                     limits = c(0, 0.45)) + #set scale size #<<
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)"
  ) +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme_barcharts +
  coord_fixed(ratio = 22, expand = TRUE, clip = "on")

taxaplot_bacter_4_nox <- taxaplot_bacter_4 +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank()
  ) +
  labs(tag = "(E)")

taxaplot_plancto_4 <- 
  data_5_ms %>% #<<
  "/"(rowSums(.)) %>% #calculate proportions
  select(names( #select names from
    which(.[1:19,] %>% #monoculture samples only
            apply(2,max) > 0) #above max value zero (select taxa with counts)
  )) %>%
  select(-which(names(.) %in% names(which(apply(., 2, max) < 0.0)))) %>% #select taxa above a percentage / remove taxa below percentage
  select(contains("Planctomycetota")) %>% #select specific taxa #<<
  rename_with(~ gsub("d__Bacteria.p__Planctomycetota.c__028H05.P.BN.P5.o__028H05.P.BN.P5.f__", "class ", .x)) %>%
  rename_with(~ gsub("d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__", "", .x)) %>%
  rename_with(~ gsub("d__Bacteria.p__Planctomycetota.c__BD7.11.o__BD7.11.f__", "class ", .x)) %>%
  rename_with(~ gsub("d__Bacteria.p__Planctomycetota.c__Phycisphaerae.o__", "", .x)) %>%
  rename_with(~ gsub("d__Bacteria.p__Planctomycetota.c__OM190.o__OM190.f__", "class ", .x)) %>%
  rename_with(~ gsub(".f__", "; ", .x)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("TB", "KF")), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_ms) %>% #<<
  mutate(SamplePlotName = fct_relevel(SamplePlotName, guni_d_order_labels_ms_dm)) %>%
  mutate(ExpGen = paste(TypeExp, Genus, sep = "_")) %>%
  ggplot(aes(x = SamplePlotName, y = abundance)) +  
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", size = 0.35) +
  scale_fill_brewer(name = "Phylum Planctomycetota (> 1%)", palette = "Paired") +
  scale_y_continuous(labels = function(x) x*100, expand = c(0, 0), limits = c(0, 0.30)) + #set scale size #<<
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)"
  ) +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme_barcharts +
  coord_fixed(ratio = 33, expand = TRUE, clip = "on")

taxaplot_plancto_4_nox <- taxaplot_plancto_4 +
  theme(#axis.title.x = element_blank(),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(tag = "(F)")

# combine plots
taxaplots_all <- 
  p_dendro_combination / taxaplot_2_nox / taxaplot_alpha_4_nox / taxaplot_gamma_4_nox / taxaplot_bacter_4_nox / taxaplot_plancto_4_nox & 
  theme(legend.justification = "left",
        plot.tag = element_text(size = 10))

# save figure and edit externally
ggsave("r_output/visualizations/mono_source_viz/taxa-plots-dendro-combo.pdf",
       plot = taxaplots_all,
       device = cairo_pdf,
       #dpi = 300,
       width = 170,
       height = 260,
       unit = "mm")
#figure editing continued in Adobe Illustator

#####----alpha div----
shannon <- read_qza("mono-source/diversity/core-metrics-results-D32260/shannon_vector.qza")
f_pd <- read_qza("mono-source/diversity/core-metrics-results-D32260/faith_pd_vector.qza")
evenness <- read_qza("mono-source/diversity/core-metrics-results-D32260/evenness_vector.qza")
observed_asv <- read_qza("mono-source/diversity/core-metrics-results-D32260/observed_features_vector.qza")

metadata_ms_alphadiv <- 
  metadata_ms %>%
  #slice_head(n = 19) %>%
  left_join(rownames_to_column(shannon$data, "SampleID")) %>%
  left_join(rownames_to_column(f_pd$data, "SampleID")) %>%
  left_join(rownames_to_column(evenness$data, "SampleID")) %>%
  left_join(rownames_to_column(observed_asv$data, "SampleID"))


#shannon_entropy <dbl>, faith_pd <dbl>, pielou_evenness <dbl>, observed_features <int>
theme_alpha <- theme(
  axis.text.x = element_text(angle = 90, vjust = 0.5, 
                             hjust = 0.9, size = 11, 
                             color = "grey30"), 
  strip.text.x = element_text(size = 11, color = "gray30", face = "bold"),
  strip.background = element_rect(color = "gray20"),
  panel.grid = element_line(colour = "gray80"),
  #panel.grid.major.x = element_line(color = "gray30", linetype = 3),
  panel.border = element_rect(color = "gray30"),
  axis.ticks.x = element_line(colour = "grey30"), # << visibility
  axis.ticks.y = element_line(colour = "grey30"),
  axis.title.x = element_text(size = 11, color = "black"),
  #axis.title.x = element_blank(), #<< visibility
  legend.text = element_markdown(size = 8, color = "grey30"),
  axis.title.y = element_text(size = 11, color = "black"),
  axis.text.y = element_text(size = 9, color = "grey30"),
  legend.title = element_text(size = 9, color = "black"),
  legend.spacing.y = unit(0.05, 'cm'),
  legend.key.size = unit(0.3, 'cm'),
  plot.tag = element_text(size = 12)
)


shannon_alpha_plot_mono <- metadata_ms_alphadiv %>%
  mutate(SamplePlotName = fct_relevel(SamplePlotName, guni_d_order_labels_ms_dm)) %>%
  mutate(Genus = fct_relevel(Genus,
                             "Achnanthes", "Amphora", "Poulinea",
                             "Diploneis", "Fallacia", "Nitzschia",
                             "Entomoneis", "Psammodictyon", "Caretta")) %>%
  filter(., TypeExp == "monoculture") %>%
  ggplot(aes(x = SamplePlotName, y = shannon_entropy)) + #<<
  geom_point(
    aes(fill = Genus, shape = TurtleID),
    alpha=1,
    size = 3,
    #shape = 21
  ) +
  geom_point(
    aes(shape = TurtleID),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Turtle ID:") +
  scale_fill_manual(name = "Genus:", 
                    values = c("#1B9E77", "#D95F02", "#7570B3", 
                               "#E7298A", "#66A61E", "#E6AB02", 
                               "#A6761D", "#666666", "#FFFFFF"),
                    labels = c(Achnanthes = "*Achnanthes*",
                               Amphora = "*Amphora*",
                               Poulinea = "*Poulinea*",
                               Diploneis = "*Diploneis*",
                               Fallacia = "*Fallacia*",
                               Nitzschia = "*Nitzschia*",
                               Entomoneis = "*Entomoneis*",
                               Psammodictyon = "*Psammodictyon*",
                               Caretta = "*Caretta* (host)")) +
  guides(shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1), 
         fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2)) +
  theme(axis.text.x = element_text(angle = 90)) +
  #facet_grid(~SourceID, scales = "free", space = "free", drop = TRUE, labeller = label_wrap_gen(width=10)) +
  labs(
    x = "Sample ID",
    y = "Shannon's entropy" #<<
  ) +
  theme_alpha

faith_alpha_plot_mono <- metadata_ms_alphadiv %>%
  mutate(SamplePlotName = fct_relevel(SamplePlotName, guni_d_order_labels_ms_dm)) %>%
  mutate(Genus = fct_relevel(Genus,
                             "Achnanthes", "Amphora", "Poulinea",
                             "Diploneis", "Fallacia", "Nitzschia",
                             "Entomoneis", "Psammodictyon", "Caretta")) %>%
  filter(., TypeExp == "monoculture") %>%
  ggplot(aes(x = SamplePlotName, y = faith_pd)) + #<<
  geom_point(
    aes(fill = Genus, shape = TurtleID),
    alpha=1,
    size = 3,
    #shape = 21
  ) +
  geom_point(
    aes(shape = TurtleID),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Turtle ID:") +
  scale_fill_manual(name = "Genus:", 
                    values = c("#1B9E77", "#D95F02", "#7570B3", 
                               "#E7298A", "#66A61E", "#E6AB02", 
                               "#A6761D", "#666666", "#FFFFFF"),
                    labels = c(Achnanthes = "*Achnanthes*",
                               Amphora = "*Amphora*",
                               Poulinea = "*Poulinea*",
                               Diploneis = "*Diploneis*",
                               Fallacia = "*Fallacia*",
                               Nitzschia = "*Nitzschia*",
                               Entomoneis = "*Entomoneis*",
                               Psammodictyon = "*Psammodictyon*",
                               Caretta = "*Caretta* (host)")) +
  guides(shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1), 
         fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2)) +
  theme(axis.text.x = element_text(angle = 90)) +
  #facet_grid(~SourceID, scales = "free", space = "free", drop = TRUE, labeller = label_wrap_gen(width=10)) +
  labs(
    x = "Sample ID",
    y = "Faith's PD" #<<
  ) +
  theme_alpha



observed_alpha_plot_mono <- metadata_ms_alphadiv %>%
  mutate(SamplePlotName = fct_relevel(SamplePlotName, guni_d_order_labels_ms_dm)) %>%
  mutate(Genus = fct_relevel(Genus,
                             "Achnanthes", "Amphora", "Poulinea",
                             "Diploneis", "Fallacia", "Nitzschia",
                             "Entomoneis", "Psammodictyon", "Caretta")) %>%
  filter(., TypeExp == "monoculture") %>%
  ggplot(aes(x = SamplePlotName, y = observed_features)) + #<<
  geom_point(
    aes(fill = Genus, shape = TurtleID),
    alpha=1,
    size = 3,
    #shape = 21
  ) +
  geom_point(
    aes(shape = TurtleID),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Turtle ID:") +
  scale_fill_manual(name = "Genus:", 
                    values = c("#1B9E77", "#D95F02", "#7570B3", 
                               "#E7298A", "#66A61E", "#E6AB02", 
                               "#A6761D", "#666666", "#FFFFFF"),
                    labels = c(Achnanthes = "*Achnanthes*",
                               Amphora = "*Amphora*",
                               Poulinea = "*Poulinea*",
                               Diploneis = "*Diploneis*",
                               Fallacia = "*Fallacia*",
                               Nitzschia = "*Nitzschia*",
                               Entomoneis = "*Entomoneis*",
                               Psammodictyon = "*Psammodictyon*",
                               Caretta = "*Caretta* (host)")) +
  guides(shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1), 
         fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2)) +
  theme(axis.text.x = element_text(angle = 90)) +
  #facet_grid(~SourceID, scales = "free", space = "free", drop = TRUE, labeller = label_wrap_gen(width=10)) +
  labs(
    x = "Sample ID",
    y = "Observed ASVs" #<<
  ) +
  theme_alpha


evenness_alpha_plot_mono <-
  metadata_ms_alphadiv %>%
  mutate(SamplePlotName = fct_relevel(SamplePlotName, guni_d_order_labels_ms_dm)) %>%
  mutate(Genus = fct_relevel(Genus,
                             "Achnanthes", "Amphora", "Poulinea",
                             "Diploneis", "Fallacia", "Nitzschia",
                             "Entomoneis", "Psammodictyon", "Caretta")) %>%
  filter(., TypeExp == "monoculture") %>%
  ggplot(aes(x = SamplePlotName, y = pielou_evenness)) + #<<
  geom_point(
    aes(fill = Genus, shape = TurtleID),
    alpha=1,
    size = 3,
    #shape = 21
  ) +
  geom_point(
    aes(shape = TurtleID),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Turtle ID:") +
  scale_fill_manual(name = "Genus:", 
                    values = c("#1B9E77", "#D95F02", "#7570B3", 
                               "#E7298A", "#66A61E", "#E6AB02", 
                               "#A6761D", "#666666", "#FFFFFF"),
                    labels = c(Achnanthes = "*Achnanthes*",
                               Amphora = "*Amphora*",
                               Poulinea = "*Poulinea*",
                               Diploneis = "*Diploneis*",
                               Fallacia = "*Fallacia*",
                               Nitzschia = "*Nitzschia*",
                               Entomoneis = "*Entomoneis*",
                               Psammodictyon = "*Psammodictyon*",
                               Caretta = "*Caretta* (host)")) +
  guides(shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1), 
         fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2)) +
  theme(axis.text.x = element_text(angle = 90)) +
  #facet_grid(~SourceID, scales = "free", space = "free", drop = TRUE, labeller = label_wrap_gen(width=10)) +
  labs(
    x = "Sample ID",
    y = "Pielou's evenness" #<<
  ) +
  theme_alpha

evenness_alpha_plot_source <- metadata_ms_alphadiv %>%
  mutate(SamplePlotName = fct_relevel(SamplePlotName, guni_d_order_labels_ms_dm)) %>%
  # mutate(Genus = fct_relevel(Genus,
  #                            "Achnanthes", "Amphora", "Poulinea",
  #                            "Diploneis", "Fallacia", "Nitzschia",
  #                            "Entomoneis", "Psammodictyon", "Caretta")) %>%
  filter(., TypeExp == "source") %>%
  ggplot(aes(x = SamplePlotName, y = pielou_evenness)) + #<<
  geom_point(
    aes(fill = Genus, shape = TurtleID),
    alpha=1,
    size = 3,
    #shape = 21
  ) +
  geom_point(
    aes(shape = TurtleID),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Turtle ID:") +
  scale_fill_manual(name = "Genus:", 
                    values = c("#FFFFFF"),
                    labels = c(Caretta = "*Caretta* (host)"))+
  guides(shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1), 
         fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2)) +
  theme(axis.text.x = element_text(angle = 90)) +
  #facet_grid(~SourceID, scales = "free", space = "free", drop = TRUE, labeller = label_wrap_gen(width=10)) +
  labs(
    x = "Sample ID",
    y = "" #<<
  ) +
  theme_alpha

faith_alpha_plot_source <- metadata_ms_alphadiv %>%
  mutate(SamplePlotName = fct_relevel(SamplePlotName, guni_d_order_labels_ms_dm)) %>%
  # mutate(Genus = fct_relevel(Genus,
  #                            "Achnanthes", "Amphora", "Poulinea",
  #                            "Diploneis", "Fallacia", "Nitzschia",
  #                            "Entomoneis", "Psammodictyon", "Caretta")) %>%
  filter(., TypeExp == "source") %>%
  ggplot(aes(x = SamplePlotName, y = faith_pd)) + #<<
  geom_point(
    aes(fill = Genus, shape = TurtleID),
    alpha=1,
    size = 3,
    #shape = 21
  ) +
  geom_point(
    aes(shape = TurtleID),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Turtle ID:") +
  scale_fill_manual(name = "Genus:", 
                    values = c("#FFFFFF"),
                    labels = c(Caretta = "*Caretta* (host)"))+
  guides(shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1), 
         fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2)) +
  theme(axis.text.x = element_text(angle = 90)) +
  #facet_grid(~SourceID, scales = "free", space = "free", drop = TRUE, labeller = label_wrap_gen(width=10)) +
  labs(
    x = "Sample ID",
    y = "" #<<
  ) +
  theme_alpha

shannon_alpha_plot_source <- metadata_ms_alphadiv %>%
  mutate(SamplePlotName = fct_relevel(SamplePlotName, guni_d_order_labels_ms_dm)) %>%
  # mutate(Genus = fct_relevel(Genus,
  #                            "Achnanthes", "Amphora", "Poulinea",
  #                            "Diploneis", "Fallacia", "Nitzschia",
  #                            "Entomoneis", "Psammodictyon", "Caretta")) %>%
  filter(., TypeExp == "source") %>%
  ggplot(aes(x = SamplePlotName, y = shannon_entropy)) + #<<
  geom_point(
    aes(fill = Genus, shape = TurtleID),
    alpha=1,
    size = 3,
    #shape = 21
  ) +
  geom_point(
    aes(shape = TurtleID),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Turtle ID:") +
  scale_fill_manual(name = "Genus:", 
                    values = c("#FFFFFF"),
                    labels = c(Caretta = "*Caretta* (host)"))+
  guides(shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1), 
         fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2)) +
  theme(axis.text.x = element_text(angle = 90)) +
  #facet_grid(~SourceID, scales = "free", space = "free", drop = TRUE, labeller = label_wrap_gen(width=10)) +
  labs(
    x = "Sample ID",
    y = "" #<<
  ) +
  theme_alpha

observed_alpha_plot_source <-
  metadata_ms_alphadiv %>%
  mutate(SamplePlotName = fct_relevel(SamplePlotName, guni_d_order_labels_ms_dm)) %>%
  # mutate(Genus = fct_relevel(Genus,
  #                            "Achnanthes", "Amphora", "Poulinea",
  #                            "Diploneis", "Fallacia", "Nitzschia",
  #                            "Entomoneis", "Psammodictyon", "Caretta")) %>%
  filter(., TypeExp == "source") %>%
  ggplot(aes(x = SamplePlotName, y = observed_features)) + #<<
  geom_point(
    aes(fill = Genus, shape = TurtleID),
    alpha=1,
    size = 3,
    #shape = 21
  ) +
  geom_point(
    aes(shape = TurtleID),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Turtle ID:") +
  scale_fill_manual(name = "Genus:", 
                    values = c("#FFFFFF"),
                    labels = c(Caretta = "*Caretta* (host)"))+
  guides(shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1), 
         fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2)) +
  theme(axis.text.x = element_text(angle = 90)) +
  #facet_grid(~SourceID, scales = "free", space = "free", drop = TRUE, labeller = label_wrap_gen(width=10)) +
  labs(
    x = "Sample ID",
    y = "" #<<
  ) +
  theme_alpha

observed_alpha_plot_mono <-
  observed_alpha_plot_mono +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

observed_alpha_plot_source <-
  observed_alpha_plot_source +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

shannon_alpha_plot_mono <-
  shannon_alpha_plot_mono +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

shannon_alpha_plot_source <-
  shannon_alpha_plot_source +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

evenness_alpha_plot_mono <-
  evenness_alpha_plot_mono +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

evenness_alpha_plot_source <-
  evenness_alpha_plot_source +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


alpha_plot <- observed_alpha_plot_mono + observed_alpha_plot_source +
  shannon_alpha_plot_mono + shannon_alpha_plot_source +
  evenness_alpha_plot_mono + evenness_alpha_plot_source +
  faith_alpha_plot_mono + faith_alpha_plot_source +
  plot_layout(widths = c(3,1), ncol = 2, guides = "collect")

ggsave(
  "r_output/visualizations/mono_source_viz/alpha-ms-shapes.pdf",
  plot = alpha_plot,
  device = cairo_pdf,
  height = 180,
  width = 168,
  units = "mm"
)


