#----------------------------------------------------------------------------	
# QIIME2 16S NGS data analyses for TurtleBIOME project (Bosak lab, University
# of Zagreb)
# By: Klara Filek
# Region: V4 region of 16S rRNA gene
# Primers: 515F and 806R (Apprill et al. 2015; Parada et al. 2016)
# Platform: Illumina MiSeq v2 (250x2 bp paired-end)
# Sequences available at European Nucleotide Archive under accessions:
# - diatom monoclonal cultures: PRJEB47668 (ERP131963)
# - source samples: PRJEB51458 (ERP136087), specifically samples ERS10917091,
# ERS10917092, ERS10917093, ERS10917103, ERS10917104, ERS10917111
# Mendeley data DOI: 10.17632/4r6568xcpw.1
#----------------------------------------------------------------------------

##----- QIIME2 environment activation in Conda ------------------------------

conda activate qiime2-2021.4

## change directory where the data is
## make directories needed for analyses

##----- Diatom monoclonal cultures NGS data (replicates) --------------------


mkdir -p mono/demux-dada2  #for demultiplexing and dada2 data
mkdir -p mono/taxonomy #for taxonomy assignment
mkdir -p mono/tree #phylogenetic trees
mkdir -p mono/filtered #data filtered based on taxonomy


## import data for diatom monoculture sequencing

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path input-data/NGS01151_Trimmed_reads_cas1.8 \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path mono/demux-dada2/demux-paired-end-monocult.qza


## summarize import stats

qiime demux summarize \
  --i-data mono/demux-dada2/demux-paired-end-monocult.qza \
  --o-visualization mono/demux-dada2/demux-paired-end-monocult.qzv


## perform DADA2 denoising (no trimming needed) 

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs mono/demux-dada2/demux-paired-end-monocult.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-n-threads 0 \
  --o-table mono/demux-dada2/table-monocult.qza \
  --o-representative-sequences mono/demux-dada2/rep-seqs-dada2-monocult.qza \
  --o-denoising-stats mono/demux-dada2/stats-dada2-monocult.qza


## summarize and tabulate to qiime visualizations

qiime metadata tabulate \
  --m-input-file mono/demux-dada2/stats-dada2-monocult.qza \
  --o-visualization mono/demux-dada2/stats-dada2-monocult.qzv

qiime feature-table summarize \
  --i-table mono/demux-dada2/table-monocult.qza \
  --o-visualization mono/demux-dada2/table-monocult.qzv \
  --m-sample-metadata-file metadata/20210510MonocultbiomeSeqMetadata.tsv

qiime feature-table tabulate-seqs \
  --i-data mono/demux-dada2/rep-seqs-dada2-monocult.qza \
  --o-visualization mono/demux-dada2/rep-seqs-dada2-monocult.qzv


## assign taxonomy with Silva silva-138-99-515-806-nb-classifier.qza 
## (MD5: e05afad0fe87542704be96ff483824d4) made with scikitlearn 0.24.1 
## downloaded from qiime2.org/2021.4/data-resources/

qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-515-806-nb-classifier.qza \
  --i-reads mono/demux-dada2/rep-seqs-dada2-monocult.qza \
  --o-classification mono/taxonomy/taxonomy.qza

qiime metadata tabulate \
  --m-input-file mono/taxonomy/taxonomy.qza \
  --o-visualization mono/taxonomy/taxonomy.qzv


## taxa bar plot prior to filtering

qiime taxa barplot \
  --i-table mono/demux-dada2/table-monocult.qza \
  --i-taxonomy mono/taxonomy/taxonomy.qza \
  --m-metadata-file metadata/20210510MonocultbiomeSeqMetadata.tsv \
  --o-visualization mono/taxonomy/taxa-bar-plots.qzv


## construct phylogeny with mafft tree

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences mono/demux-dada2/rep-seqs-dada2-monocult.qza \
  --o-alignment mono/tree/aligned-rep-seqs.qza \
  --o-masked-alignment mono/tree/masked-aligned-rep-seqs.qza \
  --o-tree mono/tree/unrooted-tree.qza \
  --o-rooted-tree mono/tree/rooted-tree.qza


## tree visualization

qiime empress community-plot \
  --i-tree mono/tree/rooted-tree.qza \
  --i-feature-table mono/demux-dada2/table-monocult.qza \
  --m-sample-metadata-file metadata/20210510MonocultbiomeSeqMetadata.tsv \
  --m-feature-metadata-file mono/taxonomy/taxonomy.qza \
  --o-visualization mono/tree/community-tree.qzv


## filtering chloroplasts, mitochondria, unassigned, eukaryota

qiime taxa filter-table \
  --i-table mono/demux-dada2/table-monocult.qza \
  --i-taxonomy mono/taxonomy/taxonomy.qza \
  --p-exclude mitochondria,chloroplast,eukaryota \
  --p-include p__ \
  --o-filtered-table mono/filtered/table-mono-fltr.qza 

qiime taxa filter-seqs \
  --i-sequences mono/demux-dada2/rep-seqs-dada2-monocult.qza \
  --i-taxonomy mono/taxonomy/taxonomy.qza \
  --p-exclude mitochondria,chloroplast,eukaryota \
  --p-include p__ \
  --o-filtered-sequences mono/filtered/rep-seqs-mono-fltr.qza

qiime feature-table summarize \
  --i-table mono/filtered/table-mono-fltr.qza \
  --o-visualization mono/filtered/table-mono-fltr.qzv \
  --m-sample-metadata-file metadata/20210510MonocultbiomeSeqMetadata.tsv

qiime feature-table tabulate-seqs \
  --i-data mono/filtered/rep-seqs-mono-fltr.qza \
  --o-visualization mono/filtered/rep-seqs-mono-fltr.qzv

## make directories for filtered data

mkdir -p mono/filtered/diversity
mkdir -p mono/filtered/taxonomy
mkdir -p mono/filtered/tree


## taxa barplot of filtered data

qiime taxa barplot \
  --i-table mono/filtered/table-mono-fltr.qza \
  --i-taxonomy mono/taxonomy/taxonomy.qza \
  --m-metadata-file metadata/20210510MonocultbiomeSeqMetadata.tsv \
  --o-visualization mono/filtered/taxonomy/taxa-bar-plots.qzv


## phylogeny mafft tree on filtered sequences

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences mono/filtered/rep-seqs-mono-fltr.qza \
  --o-alignment mono/filtered/tree/aligned-rep-seqs.qza \
  --o-masked-alignment mono/filtered/tree/masked-aligned-rep-seqs.qza \
  --o-tree mono/filtered/tree/unrooted-tree.qza \
  --o-rooted-tree mono/filtered/tree/rooted-tree.qza


## tree visualization

qiime empress community-plot \
  --i-tree mono/filtered/tree/rooted-tree.qza \
  --i-feature-table mono/filtered/table-mono-fltr.qza \
  --m-sample-metadata-file metadata/20210510MonocultbiomeSeqMetadata.tsv \
  --m-feature-metadata-file mono/taxonomy/taxonomy.qza \
  --o-visualization mono/filtered/tree/community-tree.qzv


## alpha diversity based rarefaction curves

qiime diversity alpha-rarefaction \
  --i-table mono/filtered/table-mono-fltr.qza \
  --i-phylogeny mono/filtered/tree/rooted-tree.qza \
  --p-max-depth 60000 \
  --m-metadata-file metadata/20210510MonocultbiomeSeqMetadata.tsv \
  --o-visualization mono/filtered/diversity/alpha-rarefaction.qzv


## beta diversity core metrics at 6272 depth with filtered sequences

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny mono/filtered/tree/rooted-tree.qza \
  --i-table mono/filtered/table-mono-fltr.qza \
  --p-sampling-depth 6272 \
  --m-metadata-file metadata/20210510MonocultbiomeSeqMetadata.tsv \
  --output-dir mono/filtered/diversity/core-metrics-results-D6272

## generalized unifrac (Chen et al. 2012 Bioinformatics)

mkdir mono/filtered/diversity/gen-unifrac-D6272

qiime feature-table rarefy \
  --i-table mono/filtered/table-mono-fltr.qza \
  --p-sampling-depth 6272 \
  --o-rarefied-table mono/filtered/diversity/gen-unifrac-D6272/table-D6272.qza

qiime diversity beta-phylogenetic \
  --i-table mono/filtered/diversity/gen-unifrac-D6272/table-D6272.qza \
  --i-phylogeny mono/filtered/tree/rooted-tree.qza \
  --p-metric generalized_unifrac \
  --p-threads 2 \
  --p-alpha 0.5 \
  --o-distance-matrix mono/filtered/diversity/gen-unifrac-D6272/generalized_unifrac_distance_matrix.qza

qiime diversity pcoa \
  --i-distance-matrix mono/filtered/diversity/gen-unifrac-D6272/generalized_unifrac_distance_matrix.qza \
  --o-pcoa mono/filtered/diversity/gen-unifrac-D6272/generalized_unifrac_pcoa_results.qza


## beta rarefaction tests

mkdir mono/filtered/diversity/beta-rarefaction

qiime diversity beta-rarefaction \
  --i-table mono/filtered/table-mono-fltr.qza \
  --i-phylogeny mono/filtered/tree/rooted-tree.qza\
  --p-metric generalized_unifrac \
  --p-clustering-method upgma\
  --m-metadata-file metadata/20210510MonocultbiomeSeqMetadata.tsv\
  --p-sampling-depth 6272 \
  --o-visualization mono/filtered/diversity/beta-rarefaction/betarar-gunifrac.qzv

qiime diversity beta-rarefaction \
  --i-table mono/filtered/table-mono-fltr.qza \
  --i-phylogeny mono/filtered/tree/rooted-tree.qza\
  --p-metric braycurtis \
  --p-clustering-method upgma\
  --m-metadata-file metadata/20210510MonocultbiomeSeqMetadata.tsv\
  --p-sampling-depth 6272 \
  --o-visualization mono/filtered/diversity/beta-rarefaction/betarar-bray.qzv

qiime diversity beta-rarefaction \
  --i-table mono/filtered/table-mono-fltr.qza \
  --i-phylogeny mono/filtered/tree/rooted-tree.qza\
  --p-metric aitchison \
  --p-clustering-method upgma\
  --m-metadata-file metadata/20210510MonocultbiomeSeqMetadata.tsv\
  --p-sampling-depth 6272 \
  --o-visualization mono/filtered/diversity/beta-rarefaction/betarar-aitchison.qzv

qiime diversity beta-rarefaction \
  --i-table mono/filtered/table-mono-fltr.qza \
  --i-phylogeny mono/filtered/tree/rooted-tree.qza\
  --p-metric weighted_unifrac \
  --p-clustering-method upgma\
  --m-metadata-file metadata/20210510MonocultbiomeSeqMetadata.tsv\
  --p-sampling-depth 6272 \
  --o-visualization mono/filtered/diversity/beta-rarefaction/betarar-wuni.qzv


## DEICODE (robust Aitchison) on replicate samples

mkdir mono/filtered/diversity/deicode

qiime deicode rpca \
  --i-table mono/filtered/table-mono-fltr.qza \
  --p-min-feature-count 10 \
  --p-min-sample-count 500 \
  --o-biplot mono/filtered/diversity/deicode/rAitchison-ordination.qza \
  --o-distance-matrix mono/filtered/diversity/deicode/rAitchison-distance.qza

qiime emperor biplot \
  --i-biplot mono/filtered/diversity/deicode/rAitchison-ordination.qza \
  --m-sample-metadata-file metadata/20210510MonocultbiomeSeqMetadata.tsv \
  --m-feature-metadata-file mono/taxonomy/taxonomy.qza \
  --o-visualization mono/filtered/diversity/deicode/rPCA-biplot.qzv \
  --p-number-of-features 5

mkdir mono/filtered/diversity/stats

qiime diversity beta-group-significance \
  --i-distance-matrix mono/filtered/diversity/gen-unifrac-D6272/generalized_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/20210510MonocultbiomeSeqMetadata.tsv\
  --m-metadata-column Genus \
  --p-method permanova \
  --p-pairwise True \
  --o-visualization mono/filtered/diversity/stats/guni_permanova_status.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix mono/filtered/diversity/gen-unifrac-D6272/generalized_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/20210510MonocultbiomeSeqMetadata.tsv\
  --m-metadata-column Species \
  --p-method permanova \
  --p-pairwise True \
  --o-visualization mono/filtered/diversity/stats/guni_permanova_species.qzv




##----- Diatom monoclonal cultures NGS data (grouped replicates) ------------

## Grouping samples/replicates
## diatom samples (grouping by sum is okay for these samples)

mkdir mono-grouped

qiime feature-table group \
  --i-table mono/filtered/table-mono-fltr.qza \
  --p-axis sample \
  --m-metadata-file metadata/20210510MonocultbiomeSeqMetadata.tsv \
  --m-metadata-column "Unique" \
  --p-mode sum \
  --o-grouped-table mono-grouped/table-mono.qza 

qiime feature-table summarize \
  --i-table mono-grouped/table-mono.qza \
  --o-visualization mono-grouped/table-mono.qzv \
  --m-sample-metadata-file metadata/20210510MonocultbiomeSeqGrouped.tsv


## copy mono/filtered/rep-seq-mono-fltr.qza and qzv

cp mono/filtered/rep*.qza mono-grouped

mkdir -p mono-grouped/tree

cp mono/filtered/tree/* mono-grouped/tree


## construct new tree with new metadata

qiime empress community-plot \
  --i-tree mono-grouped/tree/rooted-tree.qza \
  --i-feature-table mono-grouped/table-mono.qza \
  --m-sample-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv \
  --m-feature-metadata-file mono-grouped/taxonomy/taxonomy.qza \
  --o-visualization mono-grouped/tree/community-tree.qzv


## copy taxonomy and construct new taxa bar plots (new metadata)

mkdir -p mono-grouped/taxonomy

cp mono/taxonomy/taxonomy* mono-grouped/taxonomy

qiime taxa barplot \
  --i-table mono-grouped/table-mono.qza \
  --i-taxonomy mono-grouped/taxonomy/taxonomy.qza \
  --m-metadata-file metadata/20210510MonocultbiomeSeqGrouped.tsv \
  --o-visualization mono-grouped/taxonomy/taxa-bar-plots.qzv


## alpha diversity based rarefaction curves

mkdir -p mono-grouped/diversity

qiime diversity alpha-rarefaction \
  --i-table mono-grouped/table-mono.qza \
  --i-phylogeny mono-grouped/tree/rooted-tree.qza \
  --p-max-depth 100000 \
  --m-metadata-file metadata/20210510MonocultbiomeSeqGrouped.tsv \
  --o-visualization mono-grouped/diversity/alpha-rarefaction.qzv

## beta diversity core metrics 32260 depth with filtered and grouped sequences

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny mono-grouped/tree/rooted-tree.qza \
  --i-table mono-grouped/table-mono.qza \
  --p-sampling-depth 32260 \
  --m-metadata-file metadata/20210510MonocultbiomeSeqGrouped.tsv \
  --output-dir mono-grouped/diversity/core-metrics-results-D32260


## generalized unifrac (Chen et al. 2012 Bioinformatics)

mkdir mono-grouped/diversity/gen-unifrac-D32260

qiime feature-table rarefy \
  --i-table mono-grouped/table-mono.qza \
  --p-sampling-depth 32260 \
  --o-rarefied-table mono-grouped/diversity/gen-unifrac-D32260/table-D32260.qza

qiime diversity beta-phylogenetic \
  --i-table mono-grouped/diversity/gen-unifrac-D32260/table-D32260.qza \
  --i-phylogeny mono-grouped/tree/rooted-tree.qza \
  --p-metric generalized_unifrac \
  --p-threads 2 \
  --p-alpha 0.5 \
  --o-distance-matrix mono-grouped/diversity/gen-unifrac-D32260/generalized_unifrac_distance_matrix.qza

qiime diversity pcoa \
  --i-distance-matrix mono-grouped/diversity/gen-unifrac-D32260/generalized_unifrac_distance_matrix.qza \
  --o-pcoa mono-grouped/diversity/gen-unifrac-D32260/generalized_unifrac_pcoa_results.qza


## DEICODE (robust Aitchison) on grouped samples

mkdir mono-grouped/diversity/deicode

qiime deicode rpca \
  --i-table mono-grouped/table-mono.qza \
  --p-min-feature-count 10 \
  --p-min-sample-count 500 \
  --o-biplot mono-grouped/diversity/deicode/rAitchison-ordination.qza \
  --o-distance-matrix mono-grouped/diversity/deicode/rAitchison-distance.qza

qiime emperor biplot \
  --i-biplot mono-grouped/diversity/deicode/rAitchison-ordination.qza \
  --m-sample-metadata-file metadata/20210510MonocultbiomeSeqGrouped.tsv \
  --m-feature-metadata-file mono-grouped/taxonomy/taxonomy.qza \
  --o-visualization mono-grouped/diversity/deicode/rPCA-biplot.qzv \
  --p-number-of-features 5


## statistics

mkdir -p mono-grouped/diversity/stats


## alpha diversity group significance

## faith pd

qiime diversity alpha-group-significance \
  --i-alpha-diversity mono-grouped/diversity/core-metrics-results-D32260/faith_pd_vector.qza \
  --m-metadata-file metadata/20210510MonocultbiomeSeqGrouped.tsv \
  --o-visualization mono-grouped/diversity/stats/faith_group_sign.qzv

## shannon

qiime diversity alpha-group-significance \
  --i-alpha-diversity mono-grouped/diversity/core-metrics-results-D32260/shannon_vector.qza \
  --m-metadata-file metadata/20210510MonocultbiomeSeqGrouped.tsv \
  --o-visualization mono-grouped/diversity/stats/shannon_group_sign.qzv

## observed features

qiime diversity alpha-group-significance \
  --i-alpha-diversity mono-grouped/diversity/core-metrics-results-D32260/observed_features_vector.qza \
  --m-metadata-file metadata/20210510MonocultbiomeSeqGrouped.tsv \
  --o-visualization mono-grouped/diversity/stats/observed_featured_group_sign.qzv

## evenness

qiime diversity alpha-group-significance \
  --i-alpha-diversity mono-grouped/diversity/core-metrics-results-D32260/evenness_vector.qza \
  --m-metadata-file metadata/20210510MonocultbiomeSeqGrouped.tsv \
  --o-visualization mono-grouped/diversity/stats/evenness_group_sign.qzv


## beta diversity stats

## bray curtis
qiime diversity beta-group-significance \
  --i-distance-matrix mono-grouped/diversity/core-metrics-results-D32260/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/20210510MonocultbiomeSeqGrouped.tsv\
  --m-metadata-column TurtleID \
  --p-method permanova \
  --p-pairwise True \
  --o-visualization mono-grouped/diversity/stats/bray_curtis_permanova_turtleid.qzv

## weighted unifrac

qiime diversity beta-group-significance \
  --i-distance-matrix mono-grouped/diversity/core-metrics-results-D32260/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/20210510MonocultbiomeSeqGrouped.tsv\
  --m-metadata-column TurtleID \
  --p-method permanova \
  --p-pairwise True \
  --o-visualization mono-grouped/diversity/stats/weighted_unifrac_permanova_turtleid.qzv


##rAitchison

qiime diversity beta-group-significance \
  --i-distance-matrix mono-grouped/diversity/deicode/rAitchison-distance.qza \
  --m-metadata-file metadata/20210510MonocultbiomeSeqGrouped.tsv\
  --m-metadata-column TurtleID \
  --p-method permanova \
  --p-pairwise True \
  --o-visualization mono-grouped/diversity/stats/rAitchison_permanova_turtleid.qzv


## generalized unifrac

qiime diversity beta-group-significance \
  --i-distance-matrix mono-grouped/diversity/gen-unifrac-D32260/generalized_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/20210510MonocultbiomeSeqGrouped.tsv\
  --m-metadata-column TurtleID \
  --p-method permanova \
  --p-pairwise True \
  --o-visualization mono-grouped/diversity/stats/guni_permanova_turtleid.qzv




##----- Source samples NGS data (skin and carapace) -------------------------

mkdir -p source/data1/demux-dada2

mkdir -p source/data2/demux-dada2


## imported sequences with all non biological sequences removed


## import data1

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path input-data/NGS_env_reads_source/data1 \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path source/data1/demux-dada2/demux-paired-end-source1.qza


## import data2

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path input-data/NGS_env_reads_source/data2 \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path source/data2/demux-dada2/demux-paired-end-source2.qza


##summarize

qiime demux summarize \
  --i-data source/data1/demux-dada2/demux-paired-end-source1.qza \
  --o-visualization source/data1/demux-dada2/demux-paired-end-source1.qzv

qiime demux summarize \
  --i-data source/data2/demux-dada2/demux-paired-end-source2.qza \
  --o-visualization source/data2/demux-dada2/demux-paired-end-source2.qzv


## DADA2 denoising

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs source/data1/demux-dada2/demux-paired-end-source1.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-n-threads 0 \
  --o-table source/data1/demux-dada2/table-source1.qza \
  --o-representative-sequences source/data1/demux-dada2/rep-seqs-source1.qza \
  --o-denoising-stats source/data1/demux-dada2/denoising-stats-source1.qza

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs source/data2/demux-dada2/demux-paired-end-source2.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-n-threads 0 \
  --o-table source/data2/demux-dada2/table-source2.qza \
  --o-representative-sequences source/data2/demux-dada2/rep-seqs-source2.qza \
  --o-denoising-stats source/data2/demux-dada2/denoising-stats-source2.qza


## summarize and tabulate to qzv

qiime metadata tabulate \
  --m-input-file source/data1/demux-dada2/denoising-stats-source1.qza \
  --o-visualization source/data1/demux-dada2/denoising-stats-source1.qzv

qiime metadata tabulate \
  --m-input-file source/data2/demux-dada2/denoising-stats-source2.qza \
  --o-visualization source/data2/demux-dada2/denoising-stats-source2.qzv

qiime feature-table summarize \
  --i-table source/data1/demux-dada2/table-source1.qza \
  --o-visualization source/data1/demux-dada2/table-source1.qzv \
  --m-sample-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv

qiime feature-table summarize \
  --i-table source/data2/demux-dada2/table-source2.qza \
  --o-visualization source/data2/demux-dada2/table-source2.qzv \
  --m-sample-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv

qiime feature-table tabulate-seqs \
  --i-data source/data1/demux-dada2/rep-seqs-source1.qza \
  --o-visualization source/data1/demux-dada2/rep-seqs-source1.qzv

qiime feature-table tabulate-seqs \
  --i-data source/data2/demux-dada2/rep-seqs-source2.qza \
  --o-visualization source/data2/demux-dada2/rep-seqs-source2.qzv


## merge source samples

mkdir -p source/data-merged


## merge tables 

qiime feature-table merge \
  --i-tables source/data1/demux-dada2/table-source1.qza \
  --i-tables source/data2/demux-dada2/table-source2.qza \
  --o-merged-table source/data-merged/table-source.qza

qiime feature-table merge-seqs \
  --i-data source/data1/demux-dada2/rep-seqs-source1.qza \
  --i-data source/data2/demux-dada2/rep-seqs-source2.qza \
  --o-merged-data source/data-merged/rep-seqs-source.qza

qiime feature-table summarize \
  --i-table source/data-merged/table-source.qza \
  --o-visualization source/data-merged/table-source.qzv \
  --m-sample-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv

qiime feature-table tabulate-seqs \
  --i-data source/data-merged/rep-seqs-source.qza \
  --o-visualization source/data-merged/rep-seqs-source.qzv


## assign taxonomy

mkdir source/taxonomy

qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-515-806-nb-classifier.qza \
  --i-reads source/data-merged/rep-seqs-source.qza \
  --o-classification source/taxonomy/taxonomy-source.qza

qiime metadata tabulate \
  --m-input-file source/taxonomy/taxonomy-source.qza \
  --o-visualization source/taxonomy/taxonomy-source.qzv


## taxa bar plot unfiltered

qiime taxa barplot \
  --i-table source/data-merged/table-source.qza \
  --i-taxonomy source/taxonomy/taxonomy-source.qza \
  --m-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv \
  --o-visualization source/taxonomy/taxa-bar-plots.qzv


## phylogeny mafft tree sequences

mkdir source/tree

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences source/data-merged/rep-seqs-source.qza \
  --o-alignment source/tree/aligned-rep-seqs.qza \
  --o-masked-alignment source/tree/masked-aligned-rep-seqs.qza \
  --o-tree source/tree/unrooted-tree.qza \
  --o-rooted-tree source/tree/rooted-tree.qza


## tree visualization

qiime empress community-plot \
  --i-tree source/tree/rooted-tree.qza \
  --i-feature-table source/data-merged/table-source.qza \
  --m-sample-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv \
  --m-feature-metadata-file source/taxonomy/taxonomy-source.qza \
  --o-visualization source/tree/community-tree.qzv


## filtering

mkdir source/filtered


## filtering chloroplasts, mitochondria, unassigned, eukaryota

qiime taxa filter-table \
  --i-table source/data-merged/table-source.qza \
  --i-taxonomy source/taxonomy/taxonomy-source.qza \
  --p-exclude mitochondria,chloroplast,eukaryota,unassigned \
  --p-include p__ \
  --o-filtered-table source/filtered/table-source-fltr.qza 

qiime taxa filter-seqs \
  --i-sequences source/data-merged/rep-seqs-source.qza \
  --i-taxonomy source/taxonomy/taxonomy-source.qza \
  --p-exclude mitochondria,chloroplast,eukaryota,unassigned \
  --p-include p__ \
  --o-filtered-sequences source/filtered/rep-seqs-source-fltr.qza

qiime feature-table summarize \
  --i-table source/filtered/table-source-fltr.qza \
  --o-visualization source/filtered/table-source-fltr.qzv \
  --m-sample-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv

qiime feature-table tabulate-seqs \
  --i-data source/filtered/rep-seqs-source-fltr.qza \
  --o-visualization source/filtered/rep-seqs-source-fltr.qzv


## filtered tree and taxonomy

mkdir source/filtered/tree

mkdir source/filtered/taxonomy

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences source/filtered/rep-seqs-source-fltr.qza \
  --o-alignment source/filtered/tree/aligned-rep-seqs.qza \
  --o-masked-alignment source/filtered/tree/masked-aligned-rep-seqs.qza \
  --o-tree source/filtered/tree/unrooted-tree.qza \
  --o-rooted-tree source/filtered/tree/rooted-tree.qza


## tree visualization

qiime empress community-plot \
  --i-tree source/filtered/tree/rooted-tree.qza \
  --i-feature-table source/filtered/table-source-fltr.qza \
  --m-sample-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv \
  --m-feature-metadata-file source/taxonomy/taxonomy-source.qza \
  --o-visualization source/filtered/tree/community-tree.qzv


## taxa bar plot 

qiime taxa barplot \
  --i-table source/filtered/table-source-fltr.qza \
  --i-taxonomy source/taxonomy/taxonomy-source.qza \
  --m-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv \
  --o-visualization source/filtered/taxonomy/taxa-bar-plots.qzv


## diversity

mkdir source/filtered/diversity


## alpha rarefaction

qiime diversity alpha-rarefaction \
  --i-table source/filtered/table-source-fltr.qza \
  --i-phylogeny source/filtered/tree/rooted-tree.qza \
  --p-max-depth 100000 \
  --m-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv \
  --o-visualization source/filtered/diversity/alpha-rarefaction.qzv


## beta diversity core metrics 54000 depth with filtered and grouped sequences

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny source/filtered/tree/rooted-tree.qza \
  --i-table source/filtered/table-source-fltr.qza \
  --p-sampling-depth 54000 \
  --m-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv \
  --output-dir source/filtered/diversity/core-metrics-results-D54000




## ----- Monoculture and source data merged analyses -------------------------

mkdir mono-source


## merge tables from diatom samples (grouped) and source samples

qiime feature-table merge \
  --i-tables mono-grouped/table-mono.qza \
  --i-tables source/filtered/table-source-fltr.qza \
  --o-merged-table mono-source/table.qza

qiime feature-table merge-seqs \
  --i-data mono-grouped/rep-seqs-mono-fltr.qza \
  --i-data source/filtered/rep-seqs-source-fltr.qza \
  --o-merged-data mono-source/rep-seqs.qza

qiime feature-table summarize \
  --i-table mono-source/table.qza \
  --o-visualization mono-source/table.qzv \
  --m-sample-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv

qiime feature-table tabulate-seqs \
  --i-data mono-source/rep-seqs.qza \
  --o-visualization mono-source/rep-seqs.qzv

mkdir mono-source/taxonomy

qiime feature-table merge-taxa \
  --i-data mono-grouped/taxonomy/taxonomy.qza \
  --i-data source/taxonomy/taxonomy-source.qza \
  --o-merged-data mono-source/taxonomy/taxonomy.qza

qiime taxa barplot \
  --i-table mono-source/table.qza \
  --i-taxonomy mono-source/taxonomy/taxonomy.qza \
  --m-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv \
  --o-visualization mono-source/taxonomy/taxa-bar-plots.qzv


## construct a phylogenetic tree

mkdir mono-source/tree

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences mono-source/rep-seqs.qza \
  --o-alignment mono-source/tree/aligned-rep-seqs.qza \
  --o-masked-alignment mono-source/tree/masked-aligned-rep-seqs.qza \
  --o-tree mono-source/tree/unrooted-tree.qza \
  --o-rooted-tree mono-source/tree/rooted-tree.qza


## tree visualization

qiime empress community-plot \
  --i-tree mono-source/tree/rooted-tree.qza \
  --i-feature-table mono-source/table.qza \
  --m-sample-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv \
  --m-feature-metadata-file mono-source/taxonomy/taxonomy.qza \
  --o-visualization mono-source/tree/community-tree.qzv


## beta diversity core metrics 32260 depth with filtered and grouped sequences

mkdir mono-source/diversity

mkdir mono-source/diversity/deicode

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny mono-source/tree/rooted-tree.qza \
  --i-table mono-source/table.qza \
  --p-sampling-depth 32260 \
  --m-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv \
  --output-dir mono-source/diversity/core-metrics-results-D32260


## generalized unifrac (Chen et al. 2012 Bioinformatics)

mkdir mono-source/diversity/gen-unifrac-D32260

qiime feature-table rarefy \
  --i-table mono-source/table.qza \
  --p-sampling-depth 32260 \
  --o-rarefied-table mono-source/diversity/gen-unifrac-D32260/table-D32260.qza

qiime diversity beta-phylogenetic \
  --i-table mono-source/diversity/gen-unifrac-D32260/table-D32260.qza \
  --i-phylogeny mono-source/tree/rooted-tree.qza \
  --p-metric generalized_unifrac \
  --p-threads 2 \
  --p-alpha 0.5 \
  --o-distance-matrix mono-source/diversity/gen-unifrac-D32260/generalized_unifrac_distance_matrix.qza

qiime diversity pcoa \
  --i-distance-matrix mono-source/diversity/gen-unifrac-D32260/generalized_unifrac_distance_matrix.qza \
  --o-pcoa mono-source/diversity/gen-unifrac-D32260/generalized_unifrac_pcoa_results.qza


## DEICODE (robust Aitchison) calculations

qiime deicode rpca \
  --i-table mono-source/table.qza \
  --p-min-feature-count 10 \
  --p-min-sample-count 500 \
  --o-biplot mono-source/diversity/deicode/rAitchison-ordination.qza \
  --o-distance-matrix mono-source/diversity/deicode/rAitchison-distance.qza

qiime emperor biplot \
  --i-biplot mono-source/diversity/deicode/rAitchison-ordination.qza \
  --m-sample-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv \
  --m-feature-metadata-file mono-source/taxonomy/taxonomy.qza \
  --o-visualization mono-source/diversity/deicode/rPCA-biplot.qzv \
  --p-number-of-features 5


mkdir mono-source/diversity/beta-rarefaction

qiime diversity beta-rarefaction \
  --i-table mono-source/table.qza \
  --i-phylogeny mono-source/tree/rooted-tree.qza\
  --p-metric generalized_unifrac \
  --p-clustering-method upgma\
  --m-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv\
  --p-sampling-depth 32260 \
  --o-visualization mono-source/diversity/beta-rarefaction/betarar-guni.qzv

qiime diversity beta-rarefaction \
  --i-table mono-source/table.qza \
  --i-phylogeny mono-source/tree/rooted-tree.qza\
  --p-metric weighted_unifrac \
  --p-clustering-method upgma\
  --m-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv\
  --p-sampling-depth 32260 \
  --o-visualization mono-source/diversity/beta-rarefaction/betarar-wuni.qzv

qiime diversity beta-rarefaction \
  --i-table mono-source/table.qza \
  --i-phylogeny mono-source/tree/rooted-tree.qza\
  --p-metric aitchison \
  --p-clustering-method upgma\
  --m-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv\
  --p-sampling-depth 32260 \
  --o-visualization mono-source/diversity/beta-rarefaction/betarar-aitchison.qzv


## statistics

mkdir mono-source/diversity/stats

qiime diversity adonis \
  --i-distance-matrix mono-source/diversity/core-metrics-results-D32260/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv \
  --p-formula "Genus" \
  --p-n-jobs 2 \
  --o-visualization mono-source/diversity/stats/weighted_unifrac_adonis_genus.qzv

qiime diversity adonis \
  --i-distance-matrix mono-source/diversity/core-metrics-results-D32260/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv \
  --p-formula "Genus+TurtleID" \
  --p-n-jobs 2 \
  --o-visualization mono-source/diversity/stats/weighted_unifrac_adonis_genus-turtle.qzv

qiime diversity adonis \
  --i-distance-matrix mono-source/diversity/core-metrics-results-D32260/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv \
  --p-formula "Species+SourceID" \
  --p-n-jobs 2 \
  --o-visualization mono-source/diversity/stats/weighted_unifrac_adonis_species-source.qzv

qiime diversity adonis \
  --i-distance-matrix mono-source/diversity/core-metrics-results-D32260/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv \
  --p-formula "Species" \
  --p-n-jobs 2 \
  --o-visualization mono-source/diversity/stats/weighted_unifrac_adonis_species.qzv



##----- Merge data table, rep seq and taxonomy for monocultures and source samples -----

qiime feature-table transpose \
  --i-table mono-source/table.qza \
  --o-transposed-feature-table mono-source/transposed-table.qza

qiime metadata tabulate \
  --m-input-file mono-source/rep-seqs.qza \
  --m-input-file mono-source/taxonomy/taxonomy.qza \
  --m-input-file mono-source/transposed-table.qza \
  --o-visualization mono-source/merged-table-seq.qzv

qiime tools export \
  --input-path mono-source/merged-table-seq.qzv \
  --output-path mono-source/merged-table-seq



##----- Core features: filtering tables to desired ones ---------------------

mkdir -p core-features/collapsed_tables


## collapse tables per taxonomy

qiime taxa collapse\
  --i-table mono-source/table.qza \
  --i-taxonomy mono-source/taxonomy/taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table core-features/collapsed_tables/table-lvl2.qza

qiime taxa collapse\
  --i-table mono-source/table.qza \
  --i-taxonomy mono-source/taxonomy/taxonomy.qza \
  --p-level 3 \
  --o-collapsed-table core-features/collapsed_tables/table-lvl3.qza

qiime taxa collapse\
  --i-table mono-source/table.qza \
  --i-taxonomy mono-source/taxonomy/taxonomy.qza \
  --p-level 4 \
  --o-collapsed-table core-features/collapsed_tables/table-lvl4.qza

qiime taxa collapse\
  --i-table mono-source/table.qza \
  --i-taxonomy mono-source/taxonomy/taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table core-features/collapsed_tables/table-lvl5.qza

qiime taxa collapse\
  --i-table mono-source/table.qza \
  --i-taxonomy mono-source/taxonomy/taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table core-features/collapsed_tables/table-lvl6.qza

qiime taxa collapse\
  --i-table mono-source/table.qza \
  --i-taxonomy mono-source/taxonomy/taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table core-features/collapsed_tables/table-lvl7.qza

mkdir -p core-features/all

cp core-features/collapsed_tables/table-lvl*.qza core-features/all

cd core-features/all

for file in table-lvl*.qza; do qiime feature-table core-features --i-table ${file} --o-visualization core_all_${file}.qzv; done

rm table*.qza


## core features in mono vs source samples

mkdir -p core-features/mono-vs-source

cp metadata/20210512MonocultbiomeSeqAndSource.tsv core-features/mono-vs-source/

cp core-features/collapsed_tables/table-lvl*.qza core-features/mono-vs-source

cd core-features/mono-vs-source

for file in table-lvl*.qza; do qiime feature-table filter-samples --i-table ${file} --m-metadata-file 20210512MonocultbiomeSeqAndSource.tsv --p-where "[TypeExp]='monoculture'" --o-filtered-table mono-${file}; done

for file in table-lvl*.qza; do qiime feature-table filter-samples --i-table ${file} --m-metadata-file 20210512MonocultbiomeSeqAndSource.tsv --p-where "[TypeExp]='source'" --o-filtered-table source-${file}; done

rm table*.qza 20210512MonocultbiomeSeqAndSource.tsv

for file in mono-table-lvl*.qza; do qiime feature-table core-features --i-table ${file} --o-visualization core_${file}.qzv; done && \

for file in source-table-lvl*.qza; do qiime feature-table core-features --i-table ${file} --o-visualization core_${file}.qzv; done && 

cd ../..


## core features in epizoic vs non-epizoic

mkdir -p core-features/epizoic-vs-non

cp metadata/20210512MonocultbiomeSeqAndSource.tsv core-features/epizoic-vs-non/ 

cp core-features/collapsed_tables/table-lvl*.qza core-features/epizoic-vs-non

cd core-features/epizoic-vs-non

for file in table-lvl*.qza; do qiime feature-table filter-samples --i-table ${file} --m-metadata-file 20210512MonocultbiomeSeqAndSource.tsv --p-where "[PutativeStatus]='epizoic'" --o-filtered-table epizoic-${file}; done

for file in table-lvl*.qza; do qiime feature-table filter-samples --i-table ${file} --m-metadata-file 20210512MonocultbiomeSeqAndSource.tsv --p-where "[PutativeStatus]='non-epizoic'" --o-filtered-table non_epizoic-${file}; done

rm table*.qza 20210512MonocultbiomeSeqAndSource.tsv

for file in epizoic-table-lvl*.qza; do qiime feature-table core-features --i-table ${file} --o-visualization core_${file}.qzv; done

for file in non_epizoic-table-lvl*.qza; do qiime feature-table core-features --i-table ${file} --o-visualization core_${file}.qzv; done && cd ../..


## core features in Achnanthes genus samples

mkdir -p core-features/achnanthes

cp metadata/20210512MonocultbiomeSeqAndSource.tsv core-features/achnanthes/

cp core-features/collapsed_tables/table-lvl*.qza core-features/achnanthes

cd core-features/achnanthes

for file in table-lvl*.qza; do qiime feature-table filter-samples --i-table ${file} --m-metadata-file 20210512MonocultbiomeSeqAndSource.tsv --p-where "[Genus]='Achnanthes'" --o-filtered-table ach-${file}; done

rm table*.qza 20210512MonocultbiomeSeqAndSource.tsv

for file in ach-table-lvl*.qza; do qiime feature-table core-features --i-table ${file} --o-visualization core_${file}.qzv; done && cd ../..
