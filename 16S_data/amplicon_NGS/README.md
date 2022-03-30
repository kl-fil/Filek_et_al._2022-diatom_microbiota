# Diatom monocultures and their source samples 16S amplicon NGS data processing
Sample processing prior to sequencing and sequencing are described in detail in the **Methods** section of the Filek et al. 2022 manuscript. Primers used for sequencing were 515F and 806R (Apprill et al. 2015; Parada et al. 2016).  
Sequences were preprocessed by the sequencing company (removal of adapters, primers, non-biological seqences) while all downstream analyses were performed in QIIME 2 and R. R v4.1.1. was used for ADONIS (PERMANOVA) and data exploration and visualization.  
  
*/docs* folder contains all code and data needed for reproducing the analyses and some of the visualizations. All intermediate files produced during data analyses for the manuscript itself are available as well. QIIME 2 and R analyses can be replicated within the docs folder where all the initial files are available in folders */mono*, */mono-grouped*, */source* and */mono-source* following the *qiime_code.sh* and R scripts within the *adonis_and_dataviz.Rproj* R project. 
  
Two files were to large for github and can be found at Mendeley Data doi: 10.17632/4r6568xcpw.1; one is *demux-paired-end-monocult.qza* and the other is the SILVA classifier file *silva-138-99-515-806-nb-classifier.qza* used to assign taxonomy in this study.  
The R scripts' commands are written in a way that the output is saved in /r_output folder. R scripts contain information on the version of packages and libraries used.

### General information about QIIME 2 and the system with installed plugins is as follows:
#### System versions
Python version: 3.8.8  
QIIME 2 release: 2021.4  
QIIME 2 version: 2021.4.0  
q2cli version: 2021.4.0  
#### Installed plugins
alignment: 2021.4.0  
composition: 2021.4.0  
cutadapt: 2021.4.0  
dada2: 2021.4.0  
deblur: 2021.4.0  
deicode: 0.2.4  
demux: 2021.4.0  
diversity: 2021.4.0  
diversity-lib: 2021.4.0  
emperor: 2021.4.0  
empress: 1.2.0  
feature-classifier: 2021.4.0  
feature-table: 2021.4.0  
fragment-insertion: 2021.4.0  
gneiss: 2021.4.0  
longitudinal: 2021.4.0  
metadata: 2021.4.0  
phylogeny: 2021.4.0  
quality-control: 2021.4.0  
quality-filter: 2021.4.0  
sample-classifier: 2021.4.0  
taxa: 2021.4.0  
types: 2021.4.0  
vsearch: 2021.4.0  
