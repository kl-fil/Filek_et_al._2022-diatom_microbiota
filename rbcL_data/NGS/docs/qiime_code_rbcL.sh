#----------------------------------------------------------------------------	
# QIIME2 16S NGS data analyses for TurtleBIOME project (Bosak lab, University
# of Zagreb)
# By: Klara Filek
# Region: rbcL gene 312 barcode (Vasselon et al. 2017)
# Primers: forward-AGGTGAAYWAAAGGTTCWTAYTTAAA and reverse-CCTCTAATTTACCWACNACWG (Filek et al. 2022)
# Platform: Illumina MiSeq v2 (250x2 bp paired-end)
# Sequences available at European Nucleotide Archive under accessions:
# PRJEB51297 (ERP135907)
# Mendeley data DOI: 10.17632/4r6568xcpw.1
#----------------------------------------------------------------------------

##----- QIIME2 environment activation in Conda ------------------------------

conda activate qiime2-2021.4

## change directory where the data is
## make directories needed for analyses

##----- Source samples rbcL NGS data-----------------------------------------

## The sequnces were obtained during three different sequencing runs, therefore 
## all need to be denoised separately prior to merging
## folder input_data

## make new folders for denosised reads

mkdir -p demux_dada2/seq_round-1
mkdir -p demux_dada2/seq_round-2
mkdir -p demux_dada2/seq_round-3

##----- Diatom monoclonal cultures NGS data (replicates) --------------------


mkdir -p mono/demux-dada2  #for demultiplexing and dada2 data
mkdir -p mono/taxonomy #for taxonomy assignment
mkdir -p mono/tree #phylogenetic trees
mkdir -p mono/filtered #data filtered based on taxonomy

## Import demultiplexed samples seq_round-1, 2, and 3

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path input_data/seq_round-1 \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux_dada2/seq_round-1/demux-paired-end-1.qza

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path input_data/seq_round-2 \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux_dada2/seq_round-2/demux-paired-end-2.qza

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path input_data/seq_round-3 \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux_dada2/seq_round-3/demux-paired-end-3.qza

## summarize imported data

qiime demux summarize \
  --i-data demux_dada2/seq_round-1/demux-paired-end-1.qza \
  --o-visualization demux_dada2/seq_round-1/demux-paired-end-1.qzv

qiime demux summarize \
  --i-data demux_dada2/seq_round-2/demux-paired-end-2.qza \
  --o-visualization demux_dada2/seq_round-2/demux-paired-end-2.qzv

qiime demux summarize \
  --i-data demux_dada2/seq_round-3/demux-paired-end-3.qza \
  --o-visualization demux_dada2/seq_round-3/demux-paired-end-3.qzv


## DADA2 parameters in 2nd round 1-17 samples, no trim at 5' end for F and R sequences
## Truncate forward and reverse reads more as they could overlap with the beginning of forward reads, producing more chimera
## Minimum overlap at 30 bp
## Chimera detection "pooled"

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux_dada2/seq_round-1/demux-paired-end-1.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 160 \
  --p-trunc-len-r 160 \
  --p-min-overlap 30 \
  --p-chimera-method pooled \
  --p-n-threads 0 \
  --o-table demux_dada2/seq_round-1/table_rbcL-1.qza \
  --o-representative-sequences demux_dada2/seq_round-1/rep-seqs_rbcL-1.qza \
  --o-denoising-stats demux_dada2/seq_round-1/denoising-stats-rbcL-1.qza

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux_dada2/seq_round-2/demux-paired-end-2.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 160 \
  --p-trunc-len-r 160 \
  --p-min-overlap 30 \
  --p-chimera-method pooled \
  --p-n-threads 0 \
  --o-table demux_dada2/seq_round-2/table_rbcL-2.qza \
  --o-representative-sequences demux_dada2/seq_round-2/rep-seqs_rbcL-2.qza \
  --o-denoising-stats demux_dada2/seq_round-2/denoising-stats-rbcL-2.qza

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux_dada2/seq_round-3/demux-paired-end-3.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 160 \
  --p-trunc-len-r 160 \
  --p-min-overlap 30 \
  --p-chimera-method pooled \
  --p-n-threads 0 \
  --o-table demux_dada2/seq_round-3/table_rbcL-3.qza \
  --o-representative-sequences demux_dada2/seq_round-3/rep-seqs_rbcL-3.qza \
  --o-denoising-stats demux_dada2/seq_round-3/denoising-stats-rbcL-3.qza


## summarize DADA2 stats

qiime feature-table summarize \
  --i-table demux_dada2/seq_round-1/table_rbcL-1.qza \
  --o-visualization demux_dada2/seq_round-1/table_rbcL-1.qzv \
  --m-sample-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv

qiime feature-table summarize \
  --i-table demux_dada2/seq_round-2/table_rbcL-2.qza \
  --o-visualization demux_dada2/seq_round-2/table_rbcL-2.qzv \
  --m-sample-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv

qiime feature-table summarize \
  --i-table demux_dada2/seq_round-3/table_rbcL-3.qza \
  --o-visualization demux_dada2/seq_round-3/table_rbcL-3.qzv \
  --m-sample-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv

qiime feature-table tabulate-seqs \
  --i-data demux_dada2/seq_round-1/rep-seqs_rbcL-1.qza \
  --o-visualization demux_dada2/seq_round-1/rep-seqs_rbcL-1.qzv

qiime feature-table tabulate-seqs \
  --i-data demux_dada2/seq_round-2/rep-seqs_rbcL-2.qza \
  --o-visualization demux_dada2/seq_round-2/rep-seqs_rbcL-2.qzv

qiime feature-table tabulate-seqs \
  --i-data demux_dada2/seq_round-3/rep-seqs_rbcL-3.qza \
  --o-visualization demux_dada2/seq_round-3/rep-seqs_rbcL-3.qzv


qiime metadata tabulate \
  --m-input-file demux_dada2/seq_round-1/denoising-stats-rbcL-1.qza \
  --o-visualization demux_dada2/seq_round-1/denoising-stats-rbcL-1.qzv

qiime metadata tabulate \
  --m-input-file demux_dada2/seq_round-2/denoising-stats-rbcL-2.qza \
  --o-visualization demux_dada2/seq_round-2/denoising-stats-rbcL-2.qzv

qiime metadata tabulate \
  --m-input-file demux_dada2/seq_round-3/denoising-stats-rbcL-3.qza \
  --o-visualization demux_dada2/seq_round-3/denoising-stats-rbcL-3.qzv


## merge all runs

mkdir merged_samples

qiime feature-table merge \
  --i-tables demux_dada2/seq_round-1/table_rbcL-1.qza \
  --i-tables demux_dada2/seq_round-2/table_rbcL-2.qza \
  --i-tables demux_dada2/seq_round-3/table_rbcL-3.qza \
  --o-merged-table merged_samples/table_rbcL_merged.qza

qiime feature-table merge-seqs \
  --i-data demux_dada2/seq_round-1/rep-seqs_rbcL-1.qza \
  --i-data demux_dada2/seq_round-2/rep-seqs_rbcL-2.qza \
  --i-data demux_dada2/seq_round-3/rep-seqs_rbcL-3.qza \
  --o-merged-data merged_samples/rep-seqs_rbcL_merged.qza

qiime feature-table summarize \
  --i-table merged_samples/table_rbcL_merged.qza \
  --o-visualization merged_samples/table_rbcL_merged.qzv \
  --m-sample-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv

qiime feature-table tabulate-seqs \
  --i-data merged_samples/rep-seqs_rbcL_merged.qza \
  --o-visualization merged_samples/rep-seqs_rbcL_merged.qzv


## construct phylogenetic tree

mkdir merged_samples/tree

qiime phylogeny align-to-tree-mafft-iqtree \
  --i-sequences merged_samples/rep-seqs_rbcL_merged.qza \
  --p-n-threads auto \
  --o-alignment merged_samples/tree/aligned-rep-seqs-merged.qza \
  --o-masked-alignment merged_samples/tree/masked-aligned-rep-seqs-merged.qza \
  --o-tree merged_samples/tree/unrooted-tree-merged.qza \
  --o-rooted-tree merged_samples/tree/rooted-tree-merged.qza


## train classifier on DiatBarcode v10 database of 263bp reads (no primers, use data as is for ref seq) /DiatBarcode_v10_mot

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path DiatBarcode_v10_mot/diat_barcode_v10_263bp_mothur.fasta \
  --output-path DiatBarcode_v10_mot/diat_barcode_v10_263bp_mothur-ref-seq.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path DiatBarcode_v10_mot/diat_barcode_v10_263bp_mothur_tax.txt \
  --output-path DiatBarcode_v10_mot/ref-taxonomy.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads DiatBarcode_v10_mot/diat_barcode_v10_263bp_mothur-ref-seq.qza \
  --i-reference-taxonomy DiatBarcode_v10_mot/ref-taxonomy.qza \
  --o-classifier DiatBarcode_v10_mot/classifier-v10.qza


## import/"train" classifier and assign taxonomy

mkdir merged_samples/taxonomy

## classifier trained on MOTHUR Diat.barcode from data https://data.inrae.fr/dataset.xhtml?persistentId=doi:10.15454/V53JZV

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path diat_barcode_v10_263bp_mothur.fasta \
  --output-path diat_barcode_v10_263bp_mothur-ref-seq.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path diat_barcode_v10_263bp_mothur_tax.txt \
  --output-path ref-taxonomy.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads diat_barcode_v10_263bp_mothur-ref-seq.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier classifier-v10.qza

## assign taxonomy

qiime feature-classifier classify-sklearn \
  --i-classifier DiatBarcode_v10_mot/classifier-v10.qza \
  --i-reads merged_samples/rep-seqs_rbcL_merged.qza \
  --o-classification merged_samples/taxonomy/taxonomy.qza

qiime metadata tabulate \
  --m-input-file merged_samples/taxonomy/taxonomy.qza \
  --o-visualization merged_samples/taxonomy/taxonomy.qzv

qiime taxa barplot \
  --i-table merged_samples/table_rbcL_merged.qza \
  --i-taxonomy merged_samples/taxonomy/taxonomy.qza \
  --m-metadata-file metadata/20210512MonocultbiomeSeqAndSource.tsv \
  --o-visualization merged_samples/taxonomy/taxa-bar-plots.qzv


## Export merged reads, representative sequences and assigned taxonomy in one table

mkdir merged_samples/merged_export

qiime feature-table transpose \
  --i-table merged_samples/table_rbcL_merged.qza \
  --o-transposed-feature-table merged_samples/table_rbcL_merged-transposed.qza

qiime metadata tabulate \
  --m-int-file merged_samples/rep-seqs_rbcL_merged.qza \
  --m-input-file merged_samples/taxonomy/taxonomy.qza \
  --m-input-file merged_samples/table_rbcL_merged-transposed.qza \
  --o-visualization merged_samples/merged_export/merged-all-data.qzv

qiime tools export \
  --input-path merged_samples/merged_export/merged-all-data.qzv \
  --output-path merged_samples/merged_export/merged-all-data

