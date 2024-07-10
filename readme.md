Companion repository with helper scripts for

# OCyRS: Orthologues cyanobacterial gene-associated RNA structure motifs


This repository accompanies the main pipeline
(https://github.com/asgeissler/OCyRS-pipeline)
of the scientific publication:


"Exploring the RNA structure regulatory potential in 202 cyanobacterial genomes"  
AS  Geissler, EC Alvarez, C Anthon, NU Frigaard, J Gorodkin, and SE Seemann  
*in preparation*


## Contents

### Brenes-et-al

This folder contains script to compare by overlap the intergenic search regions
with the sRNA predicted by:

"Identification of Conserved and Potentially Regulatory Small RNAs in Heterocystous Cyanobacteria."   
Brenes-Álvarez M, Olmedo-Verd E, Vioque A, Muro-Pastor AM.    
Front Microbiol. 2016 http://journal.frontiersin.org/article/10.3389/fmicb.2016.00048

### genbank-overlap

Compare annotations from proGenomes against Rfam hits and 
GenBank, which provides the data for the
Supplementary Table S5.1 "Comparison of annotations."


### outlier trees

Plot phylogenetic trees for orthologous genes (OGs) that are in typology starkly
different from all other genes
(as determined by principle coordinate analysis in the pipeline).


### overview OG seqid

Compare the average sequence identities of orthologous genes (OGs)
to those values observed by

"A new genomic taxonomy system for the Synechococcus collective."   
Salazar VW, Tschoeke DA, Swings J, Cosenza CA, Mattoso M, Thompson CC, et al.   
Environmental Microbiology. 2020


### Public RNA-seq

The folder contains a collection of script (numerically numbered in order
of execution) that do the following steps:

0. Query NCBI database for all single-end Illumina RNA-seq datasets of any cyanobacterial species
1. Create per dataset QC pipeline job (see  https://github.com/asgeissler/RNA-Schlange)
2. Filter datasets that pass QC
3. Run per good quality datasets, the RNA-seq mapping and expression quantificaiton pipeline
   (see https://github.com/asgeissler/RNA-browser)
4. Detect expressed genes and predicted RNA structures
5. A small overview figure of species with RNA-seq data



### Reference trees

Scripts to process existing phylogenetic ideas, and to map identifiers
to a universal "NCBI taxonomy + BioProject accession + Species name"
separated by '.' symbols.

The existing trees were published in either of the following publications:

"Comparative genomics reveals insights into cyanobacterial evolution and habitat adaptation."   
Chen MY, Teng WK, Zhao L, Hu CX, Zhou YK, Han BP, et al.   
ISME J. 2021


"An Expanded Ribosomal Phylogeny of Cyanobacteria Supports a Deep Placement of Plastids."   
Moore KR.    
Frontiers in Microbiology. 2019


"ORPER: A Workflow for Constrained SSU rRNA Phylogenies."    
Cornet L, Ahn AC, Wilmotte A, Baurain D.   
Genes. 2021


### Rfam-Query

The SQL query to retrieve from the Rfam MySQL database
(https://docs.rfam.org/en/latest/database.html)
the exact number of families that have cyanobacterial sequences
in their seed alignments.


### Topology Illustration

The script generates  an illustration of the topology distance
that is used in this study. The distance metric was defined in:



"Simulation comparison of phylogeny algorithms under equal and unequal evolutionary rates."    
Kuhner, M. K. and Felsenstein, J.    
Molecular Biology and Evolution. 1994


