# Dataviz script for Filek et al. 2022 manuscript
# Processing and visualizing data related to initial sequencing results for
# monocultures (grouped replicates)
#--------------------load packages--------------------

library(tidyverse) #v1.3.1
library(ggplot2) #v3.3.5
library(qiime2R) #qiime to R import qza v0.99.6
library(patchwork) #collect multiple ggplots together v1.1.1
library(ggrepel) #repel text/labels v0.9.1
library(ggforce) #v0.3.3
library(ggtext) #v0.1.1
library(magrittr) #v2.0.1

#set global theme for ggplot
theme_set(theme_light(base_size = 11, base_family = "Lato"))

# "#<<" comment indicates the line of code that needs to be changed when using 
# a different dataset for the same code

##---- make directory for visualizations ----
dir.create("r_output/visualizations/mono_grouped_viz", recursive = TRUE)

##---- Load metadata and data-----
# load metadata in tsv format
metadata_mg <- readr::read_tsv("metadata/20210512MonocultbiomeSeqAndSource.tsv") %>%
  rename(SampleID = "#SampleID")

# load taxonomies
taxonomy_mg <- read_qza("mono-grouped/taxonomy/taxonomy.qza")

##---- Load data from qiime2 qza to R for PCoA/PCA plots-----
# load qiime pca and pcoa results.qza (in_data_16 folder)
# 16S data loading (from core-metrics-results-merged_16S-0-with-phyla-no-mitochondria-no-chloroplast-filtered-phylogeny )
# at SeqTry3 folder of 16S analyses
bray_pcoa_mg <- read_qza("mono-grouped/diversity/core-metrics-results-D32260/bray_curtis_pcoa_results.qza")
unifrac_pcoa_mg <- read_qza("mono-grouped/diversity/core-metrics-results-D32260/unweighted_unifrac_pcoa_results.qza")
w_unifrac_pcoa_mg <- read_qza("mono-grouped/diversity/core-metrics-results-D32260/weighted_unifrac_pcoa_results.qza")
jaccard_pcoa_mg <- read_qza("mono-grouped/diversity/core-metrics-results-D32260/jaccard_pcoa_results.qza")
rAitch_pca_mg <- read_qza("mono-grouped/diversity/deicode/rAitchison-ordination.qza")
g_unifrac_pcoa_mg <- read_qza("mono-grouped/diversity/gen-unifrac-D32260/generalized_unifrac_pcoa_results.qza")


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

##---- bray1_mg----
bray1_plot_mg <-
  bray_pcoa_mg$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_mg) %>%  #<<
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
  # geom_mark_ellipse(expand = 0.02, aes(color = Genus)) +
  # geom_text_repel(aes(label = SampleID), size = 1.75) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Turtle ID:") +
  scale_fill_brewer(name = "Diatom genus:", palette = "Dark2",
                    labels = c(Achnanthes = "*Achnanthes*",
                               Amphora = "*Amphora*",
                               Poulinea = "*Poulinea*",
                               Diploneis = "*Diploneis*",
                               Fallacia = "*Fallacia*",
                               Nitzschia = "*Nitzschia*",
                               Entomoneis = "*Entomoneis*",
                               Psammodictyon = "*Psammodictyon*")) +
  guides(fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical"),
         shape = guide_legend(byrow = TRUE, direction = "vertical")) +
  scale_y_continuous(breaks = c(-0.2, 0, 0.2)) +
  xlab(paste("PC1 (", round(100*bray_pcoa_mg$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*bray_pcoa_mg$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Bray-Curtis")

##---- jacc1_mg----
jacc1_plot_mg <-
  jaccard_pcoa_mg$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_mg) %>%  #<<
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
  # geom_mark_ellipse(expand = 0.02, aes(color = Genus)) +
  # geom_text_repel(aes(label = SampleID), size = 1.75) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Turtle ID:") +
  scale_fill_brewer(name = "Diatom genus:", palette = "Dark2",
                    labels = c(Achnanthes = "*Achnanthes*",
                               Amphora = "*Amphora*",
                               Poulinea = "*Poulinea*",
                               Diploneis = "*Diploneis*",
                               Fallacia = "*Fallacia*",
                               Nitzschia = "*Nitzschia*",
                               Entomoneis = "*Entomoneis*",
                               Psammodictyon = "*Psammodictyon*")) +
  guides(fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical"),
         shape = guide_legend(byrow = TRUE, direction = "vertical")) +
  xlab(paste("PC1 (", round(100*jaccard_pcoa_mg$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*jaccard_pcoa_mg$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Jaccard")

##---- wuni1_mg----
wuni1_plot_mg <-
  w_unifrac_pcoa_mg$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_mg) %>%  #<<
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
  # geom_mark_ellipse(expand = 0.02, aes(color = Genus)) +
  # geom_text_repel(aes(label = SampleID), size = 1.75) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Turtle ID:") +
  scale_fill_brewer(name = "Diatom genus:", palette = "Dark2",
                    labels = c(Achnanthes = "*Achnanthes*",
                               Amphora = "*Amphora*",
                               Poulinea = "*Poulinea*",
                               Diploneis = "*Diploneis*",
                               Fallacia = "*Fallacia*",
                               Nitzschia = "*Nitzschia*",
                               Entomoneis = "*Entomoneis*",
                               Psammodictyon = "*Psammodictyon*")) +
  guides(fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical"),
         shape = guide_legend(byrow = TRUE, direction = "vertical")) +
  xlab(paste("PC1 (", round(100*w_unifrac_pcoa_mg$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*w_unifrac_pcoa_mg$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Weighted Unifrac")

##---- uni1_mg----
uni1_plot_mg <-
  unifrac_pcoa_mg$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_mg) %>%  #<<
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
  # geom_mark_ellipse(expand = 0.02, aes(color = Genus)) +
  # geom_text_repel(aes(label = SampleID), size = 1.75) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Turtle ID:") +
  scale_fill_brewer(name = "Diatom genus:", palette = "Dark2",
                    labels = c(Achnanthes = "*Achnanthes*",
                               Amphora = "*Amphora*",
                               Poulinea = "*Poulinea*",
                               Diploneis = "*Diploneis*",
                               Fallacia = "*Fallacia*",
                               Nitzschia = "*Nitzschia*",
                               Entomoneis = "*Entomoneis*",
                               Psammodictyon = "*Psammodictyon*")) +
  guides(fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical"),
         shape = guide_legend(byrow = TRUE, direction = "vertical")) +
  xlab(paste("PC1 (", round(100*unifrac_pcoa_mg$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*unifrac_pcoa_mg$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Unweighted Unifrac")

##---- guni1_mg----
guni1_plot_mg <-
  g_unifrac_pcoa_mg$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_mg) %>%  #<<
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
  # geom_mark_ellipse(expand = 0.02, aes(color = Genus)) +
  # geom_text_repel(aes(label = SampleID), size = 1.75) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Turtle ID:") +
  scale_fill_brewer(name = "Diatom genus:", palette = "Dark2",
                    labels = c(Achnanthes = "*Achnanthes*",
                               Amphora = "*Amphora*",
                               Poulinea = "*Poulinea*",
                               Diploneis = "*Diploneis*",
                               Fallacia = "*Fallacia*",
                               Nitzschia = "*Nitzschia*",
                               Entomoneis = "*Entomoneis*",
                               Psammodictyon = "*Psammodictyon*")) +
  guides(fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical"),
         shape = guide_legend(byrow = TRUE, direction = "vertical")) +
  xlab(paste("PC1 (", round(100*g_unifrac_pcoa_mg$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*g_unifrac_pcoa_mg$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Generalized Unifrac")


##---- bray2_mg----
bray2_plot_mg <-
  bray_pcoa_mg$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_mg) %>%  #<<
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = SourceIDfig),
    alpha=0.8,
    size = 3,
    shape = 21
  ) +
  geom_point(
    size = 3,
    shape = 21,
    color = "gray30",
    fill = "transparent"
  ) +
  #geom_text_repel(aes(label = SampleID), size = 1.75) +
  # scale_fill_brewer(name = "Source sample ID:", palette = "Set1") + 
  scale_fill_manual(name = "Source sample ID:", values = c("#1F78B4", 
                                                           "#33A02C", "#B2DF8A", 
                                                           "#FF7F00", 
                                                           "#6A3D9A", "#CAB2D6")) +
  guides(fill = guide_legend(byrow = TRUE, direction = "vertical")) +
  scale_y_continuous(breaks = c(-0.2, 0, 0.2)) +
  xlab(paste("PC1 (", round(100*bray_pcoa_mg$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*bray_pcoa_mg$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #set theme
  labs(subtitle = "Bray-Curtis")

##---- jacc2_mg----
jacc2_plot_mg <-
  jaccard_pcoa_mg$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_mg) %>%  #<<
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = SourceIDfig),
    alpha=0.8,
    size = 3,
    shape = 21
  ) +
  geom_point(
    size = 3,
    shape = 21,
    color = "gray30",
    fill = "transparent"
  ) +
  #geom_text_repel(aes(label = SampleID), size = 1.75) +
  # scale_fill_brewer(name = "Source sample ID:", palette = "Set1") + 
  scale_fill_manual(name = "Source sample ID:", values = c("#1F78B4", 
                                                           "#33A02C", "#B2DF8A", 
                                                           "#FF7F00", 
                                                           "#6A3D9A", "#CAB2D6")) +
  guides(fill = guide_legend(byrow = TRUE, direction = "vertical")) +
  xlab(paste("PC1 (", round(100*jaccard_pcoa_mg$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*jaccard_pcoa_mg$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #set theme
  labs(subtitle = "Jaccard")

##---- wuni2_mg----
wuni2_plot_mg <-
  w_unifrac_pcoa_mg$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_mg) %>%  #<<
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = SourceIDfig),
    alpha=0.8,
    size = 3,
    shape = 21
  ) +
  geom_point(
    size = 3,
    shape = 21,
    color = "gray30",
    fill = "transparent"
  ) +
  #geom_text_repel(aes(label = SampleID), size = 1.75) +
  # scale_fill_brewer(name = "Source sample ID:", palette = "Set1") + 
  scale_fill_manual(name = "Source sample ID:", values = c("#1F78B4", 
                                                           "#33A02C", "#B2DF8A", 
                                                           "#FF7F00", 
                                                           "#6A3D9A", "#CAB2D6")) +
  guides(fill = guide_legend(byrow = TRUE, direction = "vertical")) +
  xlab(paste("PC1 (", round(100*w_unifrac_pcoa_mg$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*w_unifrac_pcoa_mg$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #set theme
  labs(subtitle = "Weighted Unifrac")

##---- uni2_mg----
uni2_plot_mg <-
  unifrac_pcoa_mg$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_mg) %>%  #<<
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = SourceIDfig),
    alpha=0.8,
    size = 3,
    shape = 21
  ) +
  geom_point(
    size = 3,
    shape = 21,
    color = "gray30",
    fill = "transparent"
  ) +
  #geom_text_repel(aes(label = SampleID), size = 1.75) +
  # scale_fill_brewer(name = "Source sample ID:", palette = "Set1") + 
  scale_fill_manual(name = "Source sample ID:", values = c("#1F78B4", 
                                                           "#33A02C", "#B2DF8A", 
                                                           "#FF7F00", 
                                                           "#6A3D9A", "#CAB2D6")) +
  guides(fill = guide_legend(byrow = TRUE, direction = "vertical")) +
  xlab(paste("PC1 (", round(100*unifrac_pcoa_mg$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*unifrac_pcoa_mg$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #set theme
  labs(subtitle = "Unweighted Unifrac")

##---- guni2_mg----
guni2_plot_mg <-
  g_unifrac_pcoa_mg$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_mg) %>%  #<<
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = SourceIDfig),
    alpha=0.8,
    size = 3,
    shape = 21
  ) +
  geom_point(
    size = 3,
    shape = 21,
    color = "gray30",
    fill = "transparent"
  ) +
  #geom_text_repel(aes(label = SampleID), size = 1.75) +
  # scale_fill_brewer(name = "Source sample ID:", palette = "Set1") + 
  scale_fill_manual(name = "Source sample ID:", values = c("#1F78B4", 
                                                           "#33A02C", "#B2DF8A", 
                                                           "#FF7F00", 
                                                           "#6A3D9A", "#CAB2D6")) +
  guides(fill = guide_legend(byrow = TRUE, direction = "vertical")) +
  xlab(paste("PC1 (", round(100*g_unifrac_pcoa_mg$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*g_unifrac_pcoa_mg$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #set theme
  labs(subtitle = "Generalized Unifrac")

##---- PCA plots-----
# rename column in taxonomy so we can join it to our data
taxonomy_mg <- rename(taxonomy_mg$data, FeatureID = Feature.ID)

rPCA1_plot_mg <- 
  ggplot() +
  geom_segment(data = rAitch_pca_mg$data$Species %>% #<<
                 mutate(a = sqrt(PC1^2+PC2^2+PC3^2)) %>% #calculate the distance from the origin
                 slice_max(n = 5, order_by = a) %>% #keep X furthest away points - instead of top_n -- slice_max
                 mutate(PC1=PC1*0.9, PC2=PC2*0.9) %>% #scale arrow linearly, i guess it's ok (look at emperor)
                 left_join(taxonomy_mg), #<<
               aes(x=0, xend=PC1, y=0, yend=PC2),
               size = 0.3,
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray30") +
  geom_point(
    data = rAitch_pca_mg$data$Vectors %>% #<<
      left_join(metadata_mg), #<<
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
    data = rAitch_pca_mg$data$Vectors %>% #<<
      left_join(metadata_mg), #<<
    aes(
      x = PC1, 
      y = PC2, 
      shape = TurtleID
    ),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  # geom_text_repel(data = rAitch_pca_mg$data$Vectors %>%
  #                   left_join(metadata_mg),
  #                 aes(x = PC1, y = PC2, label = SampleID), size = 1.75) +
  geom_text_repel(data = rAitch_pca_mg$data$Species %>% #<<
                    mutate(a = sqrt(PC1^2+PC2^2+PC3^2)) %>% #calculate the distance from the origin
                    slice_max(n = 5, order_by = a) %>% #keep X furthest away points - instead of top_n -- slice_max
                    mutate(PC1=PC1*0.9, PC2=PC2*0.9) %>% #scale arrow linearly, i guess it's ok (look at emperor)
                    left_join(taxonomy_mg) %>% #<<
                    mutate(Taxa = c("plain('')~italic('Alcanivorax')",
                                    "plain('')~italic('Marinobacter')",
                                    "plain('')~italic('Phycisphaera')",
                                    "plain('')~italic('Neptuniibacter')",
                                    "plain('')~italic('Alteromonas')")),
                  aes(x = PC1, y = PC2, 
                      label = Taxa),
                  point.padding = 2,
                  alpha = 0.8,
                  colour = "black",
                  size = 3,
                  parse = TRUE) +
  scale_shape_manual(name = "Turtle ID:", values = c(21, 22, 23, 24)) +
  scale_fill_brewer(name = "Diatom genus:", palette = "Dark2",
                    labels = c(Achnanthes = "*Achnanthes*",
                               Amphora = "*Amphora*",
                               Poulinea = "*Poulinea*",
                               Diploneis = "*Diploneis*",
                               Fallacia = "*Fallacia*",
                               Nitzschia = "*Nitzschia*",
                               Entomoneis = "*Entomoneis*",
                               Psammodictyon = "*Psammodictyon*")) +
  guides(shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1), 
         fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2)) +
  xlab(paste("PC1 (",round(100*rAitch_pca_mg$data$ProportionExplained[1],2),"%)")) + #<<
  ylab(paste("PC2 (",round(100*rAitch_pca_mg$data$ProportionExplained[2],2),"%)")) + #<<
  theme_pcoa1 +
  labs(subtitle = "Robust Aitchison")

rPCA2_plot_mg <-
  rAitch_pca_mg$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_mg) %>%  #<<
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = SourceIDfig),
    alpha=0.8,
    size = 3,
    shape = 21
  ) +
  geom_point(
    size = 3,
    shape = 21,
    color = "gray30",
    fill = "transparent"
  ) +
  #geom_text_repel(aes(label = SampleID), size = 1.75) +
  #scale_fill_brewer(name = "Source sample ID:", palette = "Set1") + 
  scale_fill_manual(name = "Source sample ID:", values = c("#1F78B4", 
                                                           "#33A02C", "#B2DF8A", 
                                                           "#FF7F00", 
                                                           "#6A3D9A", "#CAB2D6")) +
  guides(fill = guide_legend(byrow = TRUE, direction = "vertical")) +
  xlab(paste("PC1 (", round(100*rAitch_pca_mg$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*rAitch_pca_mg$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 +
  labs(subtitle = "Robust Aitchison")

# plots

diss_manuscript_plot_mg <- (rPCA1_plot_mg + guni1_plot_mg) /
  (rPCA2_plot_mg + guni2_plot_mg) +
  plot_annotation(tag_levels = 'A', tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect")

diss_genus_plot_mg <- (bray1_plot_mg + wuni1_plot_mg) / (jacc1_plot_mg + uni1_plot_mg) +
  plot_annotation(tag_levels = 'A', tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect")

diss_source_plot_mg <- (bray2_plot_mg + wuni2_plot_mg) / (jacc2_plot_mg + uni2_plot_mg) +
  plot_annotation(tag_levels = 'A', tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect")

ggsave("r_output/visualizations/mono_grouped_viz/manuscript_fig_dissim_guni_rpca.pdf",
       plot = diss_manuscript_plot_mg,
       #dpi = ,
       device = cairo_pdf,
       width = 168,
       height = 140, 
       units = c("mm"))

ggsave("r_output/visualizations/mono_grouped_viz/suppl_fig_dissim_rest_genus.pdf",
       plot = diss_genus_plot_mg,
       #dpi = ,
       device = cairo_pdf,
       width = 168,
       height = 140, 
       units = c("mm"))

ggsave("r_output/visualizations/mono_grouped_viz/suppl_fig_dissim_rest_source.pdf",
       plot = diss_source_plot_mg,
       #dpi = ,
       device = cairo_pdf,
       width = 168,
       height = 140, 
       units = c("mm"))


##---- saving pcoa/pca plots -----
# add margins to the right of the plots (for 1st column in multipanel plots)

bray1_plot_mg_a <- bray1_plot_mg +
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))

jacc1_plot_mg_a <- jacc1_plot_mg +
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))

wuni1_plot_mg_a <- wuni1_plot_mg +
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))

uni1_plot_mg_a <- uni1_plot_mg +
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))

guni1_plot_mg_a <- guni1_plot_mg +
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))

rPCA1_plot_mg_a <- rPCA1_plot_mg +
  theme(plot.margin = unit(c(0,0.8,0,0), "cm"))

all_ord_plots2_mg_a <- (bray1_plot_mg_a + bray2_plot_mg + jacc1_plot_mg_a + jacc2_plot_mg) /
  (wuni1_plot_mg_a + wuni2_plot_mg + uni1_plot_mg_a + uni2_plot_mg) /
  (guni1_plot_mg_a + guni2_plot_mg + rPCA1_plot_mg_a + rPCA2_plot_mg) +
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        plot.tag.position = c(0.015,0.995),
        plot.tag = element_text(size = 10))

ggsave("r_output/visualizations/mono_grouped_viz/all-250-600-test.pdf",
       plot = all_ord_plots2_mg_a,
       device = cairo_pdf, #cairo_pdf, png or tiff
       #dpi = 350, #uncomment for png, tiff/ comment for pdf
       width = 250, 
       height = 600, #240 long 150 short
       units = c("mm"))

