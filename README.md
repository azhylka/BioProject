BioProject
==========

Contains implemented methods and algorithms used in the research on analysing micobacterium tuberculosis drug resistance.

### Important files:

##### filtering.R
Contains filtering methods:
* Minor allele frequency test

##### hmm.R
Contains adoption of the algorithm that uses Hidden Markov Models (HMM) for testing FDR and described by Zhi We in  "Hidden Markov Models for Controlling False Discovery Rate in Genome-Wide Association Analysis" [article](http://www.springerprotocols.com/Abstract/doi/10.1007/978-1-61779-400-1_22). It also makes use of PLIS package implemented by Zhi Wei himself in R.

##### distances_matrix.py
Used to build matrix of Hamming distances between all the provided genomes.

##### clustering.R
Provides clustering of the SNPs set based on the provided set of significant SNPs. That means that each cluster contains at least one significant SNP and all those SNPs that have occur wherever the target significant SNP does. This approach is currently being under consideration of reconstruction because of its results do not provide the desirable efficiency.

##### distribution.R
Calculates probabilities of the SNP being the cause of the resistance based on the set of clusters. Consequently, this module is dependent on clustering.R module so that its results currently should not be judged nor considered complete or adequate.

##### mutation_distances.Rmd
R Markdown file that was used in order to analyse possible dependencies between the distances between the genomes and their phenotype.

### Util files

##### setup.R
Contains code used to setup the working environment, e.g. working directory, common imports.

##### converter.ipynb
Used to convert some genuine files to more convenient form, like apr15.snps.matrix, which still will need to be transposed.
