core-genome-alignment
=====================

Scripts to make a concatenated alignment of core genes from bacterial genomes.

###RunProkka.ipy
This script runs the annotations program Prokka on a directory of fasta format genomes. Each genome will hava a directory as output. Additionally, a directory containing protein sequences in fasta format for each genome will be created as well as a directory containing the corresponding nucleotide sequences.

Requirements: IPython, Prokka (http://www.vicbioinformatics.com/software.prokka.shtml)

Current Versions: IPython v 0.13.1, Prokka v 1.7

Usage: RunProkka.ipy [fasta directory]

###RunOrthoMCL.ipy

This script runs the program OrthoMCL to put proteins from annotated genomes
into orthologous groups. A sql database must be created and empty for this
program to run correctly. Details of the database should be provided in the   
OrthoMCL config file. Additionally, the directory where the annotated proteins
are located should be provided to this script. Do not use relative paths.

Requirements: IPython, OrthoMCL (http://orthomcl.org/orthomcl/)

Current Versions: IPython v 0.13.1, OrthoMCL v 2.0.9

Usage: RunOrthoMCL.ipy [protein fasta directory] [OrthoMCL config file]

###FilterOrthoMCLGroups.py

This script filters protein groups output from OrthoMCL and returns those that contain only one protein per genome and have all genomes represented (the core proteins).

Usage: FilterOrthoMCLGroups.py -g [input groups file] -n [number of genomes]

###MakeCoreGenomeAlignment.ipy

This script makes a core genome alignment from proteins in a groups file in the format output by OrthoMCL. FilterOrthoMCLGroups.py should be used first. The script should be provided with fasta files of DNA sequences for proteins located in one directory.

Requirements: IPython, Biopython, MAFFT, trimAl

Current Versions: IPython v 0.13.1, Biopython v 1.63, MAFFT v 7.130b, trimAl v 1.3

Usage: MakeCoreGenomeAlignment.ipy [input groups file] [directory with DNA sequences]

###GeneMLTrees.ipy

This script produces maximum likelihood trees for a group of gene alignments
stored in a directory. The script is designed to use output from
MakeCoreGenomeAlignment.ipy, and it looks for alignments with file names
ending in _trim.fasta. The script uses RAxML for phylogenetic inference.
The scripts should run within directory containing gene alignments

Requirements: IPython, RAxML (https://github.com/stamatak/standard-RAxML)

Current Versions: IPython 0.13.1, RAxML 8.0.6

Usage: GeneMLTrees.ipy

###mdsRFdistance.R

This script produces a plot of multidimensional scaling of Robinson-Foulds
distances produced by RAxML.

Requirements: ggplot2, reshape2

Usage: mdsRFdistance.R [input file]

###eggNOGqueryFile.py

This script takes the first protein from each file in a directory of unaligned
fasta files and changes the protein ID to the name of the file. All protein  
sequences are written to an output file that should be used as the query for
eggNOGblast.py.

Requirements: Biopython

Current Versions: Python 2.7.3, Biopython v 1.63

Usage: eggNOGqueryFile.py [directory of FASTA protein sequences]

###eggNOGblast.py

This script takes a file with fasta genes and uses NCBI BLAST+ to look for
matches in the eggNOG BLAST database. It also uses eggNOG to assign functional
categories to each fasta. Inputs are the directory where eggnog files are
found and an input fasta.

Requirements: Biopython, NCBI BLAST+ (http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download), eggNOG files (http://eggnog.embl.de/version_4.0.beta/downloads.v4.html), BLAST database created from eggNOG proteins

Current Versions: Python 2.7.3, Biopython 1.63, BLAST+ 2.2.18, eggNOG 4

Usage: eggNOGblast.py -f [input fasta file] -e [directory with eggNOG files]

###geneContentTree.py

This script uses the orthologous groups output by OrthoMCL to create a    
distance matrix between strains. A neighbor joining tree is created from this
distance matrix.

Requirements: Biopython

Current Versions: Python 2.7.3, Biopython 1.64 (unreleased from https://github.com/biopython/biopython)

Usage: geneContentTree.py -g [OrthoMCL groups file]

###geneContentMatrix.py

This script uses the orthologous groups output by OrthoMCL to create a
matrix of gene content for each genome in the analysis. Each row in the
matrix is an orthologous group, and each column is a genome. The numbers
correspond to the number of members of the group contained within each
genome.

Current Versions: Python 2.7.3

Usage: geneContentMatrix.py -g [OrthoMCL groups file]

###compareCoreGenomes.py

This script used the orthologous groups output by OrthoMCL to create files
describing the core genome and genes unique to the core genome of particular
clades.

Format of clade file:

StrainA Clade1

StrainB Clade1

StrainC Clade2

Current Versions: Python 2.7.3

Usage: compareCoreGenomes.py -g [OrthoMCL groups file] -c [clades file]
