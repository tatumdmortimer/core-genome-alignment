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
