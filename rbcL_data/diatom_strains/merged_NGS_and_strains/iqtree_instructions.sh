#Iqtree2 installation 
#tutorial: http://www.iqtree.org/doc/Tutorial

$ conda install -c bioconda iqtree

# run iqtree2 on fasta alignment
# use -nt AUTO for automatic determining of how many threads the program uses
# use model finder to determine the best substitution model automatically (default)

# cite: S. Kalyaanamoorthy, B.Q. Minh, T.K.F. Wong, A. von Haeseler, and L.S. Jermiin
# (2017) ModelFinder: fast model selection for accurate phylogenetic estimates. Nat. 
# Methods, 14:587–589. DOI: 10.1038/nmeth.4285

#use ultrafast bootstrap -bb 1000


$ iqtree -s ngs_mono_merged-alignment.fas -m MFP -bb 1000 -nt AUTO

