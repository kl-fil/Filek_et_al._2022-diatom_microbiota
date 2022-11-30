# ADONIS tests script for Filek et al. 2022 manuscript

library(qiime2R) #qiime to R import qza v0.99.6
library(tidyverse) #v1.3.1
library(vegan) #2.5.7
library(pairwiseAdonis) #v0.4
library(magrittr) #v2.0.1

set.seed(19121216)

dir.create(path = "r_output/adonis_stats", recursive = TRUE)

# import metadata
metadata_ms <- read_tsv("metadata/20210512MonocultbiomeSeqAndSource.tsv") %>%
  rename(SampleID = "#SampleID")

metadata_mg <- read_tsv("metadata/20210510MonocultbiomeSeqGrouped.tsv") %>%
  rename(SampleID = "#SampleID")

metadata_m <- read_tsv("metadata/20210510MonocultbiomeSeqMetadata.tsv") %>%
  rename(SampleID = "#SampleID")

# see variable names
colnames(metadata_ms)

## import distances and extract data
## bray-curtis----
bray_distances_m <- read_qza("mono/filtered/diversity/core-metrics-results-D6272/bray_curtis_distance_matrix.qza")
bray_distances_mg <- read_qza("mono-grouped/diversity/core-metrics-results-D32260/bray_curtis_distance_matrix.qza")
bray_distances_ms <- read_qza("mono-source/diversity/core-metrics-results-D32260/bray_curtis_distance_matrix.qza")

bray_dist_m <- bray_distances_m$data
bray_dist_mg <- bray_distances_mg$data
bray_dist_ms <- bray_distances_ms$data

## jaccard----
jacc_distances_m <- read_qza("mono/filtered/diversity/core-metrics-results-D6272/jaccard_distance_matrix.qza")
jacc_distances_mg <- read_qza("mono-grouped/diversity/core-metrics-results-D32260/jaccard_distance_matrix.qza")
jacc_distances_ms <- read_qza("mono-source/diversity/core-metrics-results-D32260/jaccard_distance_matrix.qza")

jacc_dist_m <- jacc_distances_m$data
jacc_dist_mg <- jacc_distances_mg$data
jacc_dist_ms <- jacc_distances_ms$data

## generalized unifrac----
gen_unifrac_distances_ms <- read_qza("mono-source/diversity/gen-unifrac-D32260/generalized_unifrac_distance_matrix.qza")
gen_unifrac_distances_mg <- read_qza("mono-grouped/diversity/gen-unifrac-D32260/generalized_unifrac_distance_matrix.qza")
gen_unifrac_distances_m <- read_qza("mono/filtered/diversity/gen-unifrac-D6272/generalized_unifrac_distance_matrix.qza")

guni_dist_ms <- gen_unifrac_distances_ms$data
guni_dist_mg <- gen_unifrac_distances_mg$data
guni_dist_m <- gen_unifrac_distances_m$data

## unweighted unifrac----
uni_unifrac_distances_ms <- read_qza("mono-source/diversity/core-metrics-results-D32260/unweighted_unifrac_distance_matrix.qza")
uni_unifrac_distances_mg <- read_qza("mono-grouped/diversity/core-metrics-results-D32260/unweighted_unifrac_distance_matrix.qza")
uni_unifrac_distances_m <- read_qza("mono/filtered/diversity/core-metrics-results-D6272/unweighted_unifrac_distance_matrix.qza")

uni_dist_ms <- uni_unifrac_distances_ms$data
uni_dist_mg <- uni_unifrac_distances_mg$data
uni_dist_m <- uni_unifrac_distances_m$data

## weighted unifrac----
wuni_unifrac_distances_ms <- read_qza("mono-source/diversity/core-metrics-results-D32260/weighted_unifrac_distance_matrix.qza")
wuni_unifrac_distances_mg <- read_qza("mono-grouped/diversity/core-metrics-results-D32260/weighted_unifrac_distance_matrix.qza")
wuni_unifrac_distances_m <- read_qza("mono/filtered/diversity/core-metrics-results-D6272/weighted_unifrac_distance_matrix.qza")

wuni_dist_ms <- wuni_unifrac_distances_ms$data
wuni_dist_mg <- wuni_unifrac_distances_mg$data
wuni_dist_m <- wuni_unifrac_distances_m$data

## robust Aitchison----
raitch_unifrac_distances_ms <- read_qza("mono-source/diversity/deicode/rAitchison-distance.qza")
raitch_unifrac_distances_mg <- read_qza("mono-grouped/diversity/deicode/rAitchison-distance.qza")
raitch_unifrac_distances_m <- read_qza("mono/filtered/diversity/deicode/rAitchison-distance.qza")

raitch_dist_ms <- raitch_unifrac_distances_ms$data
raitch_dist_mg <- raitch_unifrac_distances_mg$data
raitch_dist_m <- raitch_unifrac_distances_m$data

## adonis runs and export ----
# group replicate monoculture samples by Source sample
adonis_rep_mono_source <- list()
adonis_rep_mono_source[["Bray"]] <- as.data.frame(adonis(bray_dist_m~SourceID, data = metadata_m, permutations = 9999)[[1]])
adonis_rep_mono_source[["Jacc"]] <- as.data.frame(adonis(jacc_dist_m~SourceID, data = metadata_m, permutations = 9999)[[1]])
adonis_rep_mono_source[["gUni"]] <- as.data.frame(adonis(guni_dist_m~SourceID, data = metadata_m, permutations = 9999)[[1]])
adonis_rep_mono_source[["uUni"]] <- as.data.frame(adonis(uni_dist_m~SourceID, data = metadata_m, permutations = 9999)[[1]])
adonis_rep_mono_source[["wUni"]] <- as.data.frame(adonis(wuni_dist_m~SourceID, data = metadata_m, permutations = 9999)[[1]])
adonis_rep_mono_source[["rAitch"]] <- as.data.frame(adonis(raitch_dist_m~SourceID, data = metadata_m, permutations = 9999)[[1]])

adonis_rep_mono_source <- do.call(rbind, adonis_rep_mono_source)

write.csv(adonis_rep_mono_source, "r_output/adonis_stats/adonis_monocult_replicates_by_source_sample.csv")


# group replicate monoculture samples TurtleID
adonis_rep_mono_turtleid <- list()
adonis_rep_mono_turtleid[["Bray"]] <- as.data.frame(adonis(bray_dist_m~TurtleID, data = metadata_m, permutations = 9999)[[1]])
adonis_rep_mono_turtleid[["Jacc"]] <- as.data.frame(adonis(jacc_dist_m~TurtleID, data = metadata_m, permutations = 9999)[[1]])
adonis_rep_mono_turtleid[["gUni"]] <- as.data.frame(adonis(guni_dist_m~TurtleID, data = metadata_m, permutations = 9999)[[1]])
adonis_rep_mono_turtleid[["uUni"]] <- as.data.frame(adonis(uni_dist_m~TurtleID, data = metadata_m, permutations = 9999)[[1]])
adonis_rep_mono_turtleid[["wUni"]] <- as.data.frame(adonis(wuni_dist_m~TurtleID, data = metadata_m, permutations = 9999)[[1]])
adonis_rep_mono_turtleid[["rAitch"]] <- as.data.frame(adonis(raitch_dist_m~TurtleID, data = metadata_m, permutations = 9999)[[1]])

adonis_rep_mono_turtleid <- do.call(rbind, adonis_rep_mono_turtleid)

write.csv(adonis_rep_mono_turtleid, "r_output/adonis_stats/adonis_monocult_replicates_by_turtleid.csv")

# group monoculture samples by TurtleID
adonis_mono_turtleid <- list()
adonis_mono_turtleid[["Bray"]] <- as.data.frame(adonis(bray_dist_mg~TurtleID, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_turtleid[["Jacc"]] <- as.data.frame(adonis(jacc_dist_mg~TurtleID, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_turtleid[["gUni"]] <- as.data.frame(adonis(guni_dist_mg~TurtleID, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_turtleid[["uUni"]] <- as.data.frame(adonis(uni_dist_mg~TurtleID, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_turtleid[["wUni"]] <- as.data.frame(adonis(wuni_dist_mg~TurtleID, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_turtleid[["rAitch"]] <- as.data.frame(adonis(raitch_dist_mg~TurtleID, data = metadata_mg, permutations = 9999)[[1]])

adonis_mono_turtleid <- do.call(rbind, adonis_mono_turtleid)

write.csv(adonis_mono_turtleid, "r_output/adonis_stats/adonis_monocultures_by_turtleid.csv")

adonis_pair_mono_turtleid <- list()
adonis_pair_mono_turtleid[["Bray"]] <- pairwise.adonis(bray_dist_mg, metadata_mg$TurtleID)
adonis_pair_mono_turtleid[["Jacc"]] <- pairwise.adonis(jacc_dist_mg, metadata_mg$TurtleID)
adonis_pair_mono_turtleid[["gUni"]] <- pairwise.adonis(guni_dist_mg, metadata_mg$TurtleID)
adonis_pair_mono_turtleid[["uUni"]] <- pairwise.adonis(uni_dist_mg, metadata_mg$TurtleID)
adonis_pair_mono_turtleid[["wUni"]] <- pairwise.adonis(wuni_dist_mg, metadata_mg$TurtleID)
adonis_pair_mono_turtleid[["rAitch"]] <- pairwise.adonis(raitch_dist_mg, metadata_mg$TurtleID)

adonis_pair_mono_turtleid <- do.call(rbind, adonis_pair_mono_turtleid)

write.csv(adonis_pair_mono_turtleid, "r_output/adonis_stats/adonis_monocultures_by_turtleid_pairwise.csv")

# group monoculture samples by SourceID
adonis_mono_sourceid <- list()
adonis_mono_sourceid[["Bray"]] <- as.data.frame(adonis(bray_dist_mg~SourceID, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_sourceid[["Jacc"]] <- as.data.frame(adonis(jacc_dist_mg~SourceID, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_sourceid[["gUni"]] <- as.data.frame(adonis(guni_dist_mg~SourceID, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_sourceid[["uUni"]] <- as.data.frame(adonis(uni_dist_mg~SourceID, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_sourceid[["wUni"]] <- as.data.frame(adonis(wuni_dist_mg~SourceID, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_sourceid[["rAitch"]] <- as.data.frame(adonis(raitch_dist_mg~SourceID, data = metadata_mg, permutations = 9999)[[1]])

adonis_mono_sourceid <- do.call(rbind, adonis_mono_sourceid)

write.csv(adonis_mono_sourceid, "r_output/adonis_stats/adonis_monocultures_by_sourceid.csv")

pairwise.adonis(raitch_dist_mg, metadata_mg$Genus)

# group monoculture samples by putative status (epizoic vs non-epizoic)
adonis_mono_putativestatus <- list()
adonis_mono_putativestatus[["Bray"]] <- as.data.frame(adonis(bray_dist_mg~PutativeStatus, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_putativestatus[["Jacc"]] <- as.data.frame(adonis(jacc_dist_mg~PutativeStatus, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_putativestatus[["gUni"]] <- as.data.frame(adonis(guni_dist_mg~PutativeStatus, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_putativestatus[["uUni"]] <- as.data.frame(adonis(uni_dist_mg~PutativeStatus, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_putativestatus[["wUni"]] <- as.data.frame(adonis(wuni_dist_mg~PutativeStatus, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_putativestatus[["rAitch"]] <- as.data.frame(adonis(raitch_dist_mg~PutativeStatus, data = metadata_mg, permutations = 9999)[[1]])

adonis_mono_putativestatus <- do.call(rbind, adonis_mono_putativestatus)

write.csv(adonis_mono_putativestatus, "r_output/adonis_stats/adonis_monocultures_by_putative_status.csv")


# group monoculture samples by TurtleID and Genus
adonis_mono_turtle_genus <- list()
adonis_mono_turtle_genus[["Bray"]] <- as.data.frame(adonis(bray_dist_mg~TurtleID+Genus, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_turtle_genus[["Jacc"]] <- as.data.frame(adonis(jacc_dist_mg~TurtleID+Genus, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_turtle_genus[["gUni"]] <- as.data.frame(adonis(guni_dist_mg~TurtleID+Genus, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_turtle_genus[["uUni"]] <- as.data.frame(adonis(uni_dist_mg~TurtleID+Genus, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_turtle_genus[["wUni"]] <- as.data.frame(adonis(wuni_dist_mg~TurtleID+Genus, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_turtle_genus[["rAitch"]] <- as.data.frame(adonis(raitch_dist_mg~TurtleID+Genus, data = metadata_mg, permutations = 9999)[[1]])

adonis_mono_turtle_genus <- do.call(rbind, adonis_mono_turtle_genus)

write.csv(adonis_mono_turtle_genus, "r_output/adonis_stats/adonis_monocultures_by_turtle_and_genus.csv")

# group monoculture samples by SourceID and Genus
adonis_mono_source_genus <- list()
adonis_mono_source_genus[["Bray"]] <- as.data.frame(adonis(bray_dist_mg~SourceID+Genus, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_source_genus[["Jacc"]] <- as.data.frame(adonis(jacc_dist_mg~SourceID+Genus, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_source_genus[["gUni"]] <- as.data.frame(adonis(guni_dist_mg~SourceID+Genus, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_source_genus[["uUni"]] <- as.data.frame(adonis(uni_dist_mg~SourceID+Genus, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_source_genus[["wUni"]] <- as.data.frame(adonis(wuni_dist_mg~SourceID+Genus, data = metadata_mg, permutations = 9999)[[1]])
adonis_mono_source_genus[["rAitch"]] <- as.data.frame(adonis(raitch_dist_mg~SourceID+Genus, data = metadata_mg, permutations = 9999)[[1]])

adonis_mono_source_genus <- do.call(rbind, adonis_mono_source_genus)

write.csv(adonis_mono_source_genus, "r_output/adonis_stats/adonis_monocultures_by_source_and_genus.csv")

# group all samples by TypeExp (monoculture vs source)
adonis_all_typeexp <- list()
adonis_all_typeexp[["Bray"]] <- as.data.frame(adonis(bray_dist_ms~TypeExp, data = metadata_ms, permutations = 9999)[[1]])
adonis_all_typeexp[["Jacc"]] <- as.data.frame(adonis(jacc_dist_ms~TypeExp, data = metadata_ms, permutations = 9999)[[1]])
adonis_all_typeexp[["gUni"]] <- as.data.frame(adonis(guni_dist_ms~TypeExp, data = metadata_ms, permutations = 9999)[[1]])
adonis_all_typeexp[["uUni"]] <- as.data.frame(adonis(uni_dist_ms~TypeExp, data = metadata_ms, permutations = 9999)[[1]])
adonis_all_typeexp[["wUni"]] <- as.data.frame(adonis(wuni_dist_ms~TypeExp, data = metadata_ms, permutations = 9999)[[1]])
adonis_all_typeexp[["rAitch"]] <- as.data.frame(adonis(raitch_dist_ms~TypeExp, data = metadata_ms, permutations = 9999)[[1]])

adonis_all_typeexp <- do.call(rbind, adonis_all_typeexp)

write.csv(adonis_all_typeexp, "r_output/adonis_stats/adonis_all_by_type_exp.csv")
