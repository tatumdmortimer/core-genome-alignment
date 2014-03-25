core-genome-alignment
=====================

Scripts to make a concatenated alignment of core genes from bacterial genomes.

##RunProkka.ipy
This script runs the annotations program Prokka on a directory of fasta format genomes. Each genome will hava a directory as output. Additionally, a directory containing protein sequences in fasta format for each genome will be created as well as a directory containing the corresponding nucleotide sequences.

Requirements: IPython, Prokka (http://www.vicbioinformatics.com/software.prokka.shtml)

Current Versions: IPython v 0.13.1, Prokka v 1.7

Usage: RunProkka.ipy [fasta directory]

##RunOrthoMCL.ipy

This script runs the program OrthoMCL to put proteins from annotated genomes
into orthologous groups. A sql database must be created and empty for this
program to run correctly. Details of the database should be provided in the   
OrthoMCL config file. Additionally, the directory where the annotated proteins
are located should be provided to this script. Do not use relative paths.

Requirements: IPython, OrthoMCL (http://orthomcl.org/orthomcl/)

Current Versions: IPython v 0.13.1, OrthoMCL v 2.0.9

Usage: RunOrthoMCL.ipy [protein fasta directory] [OrthoMCL config file]
