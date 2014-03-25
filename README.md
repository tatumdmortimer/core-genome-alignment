core-genome-alignment
=====================

Scripts to make a concatenated alignment of core genes from bacterial genomes.

##RunProkka.ipy
This script runs the annotations program Prokka on a directory of fasta format genomes. Each genome will hava a directory as output. Additionally, a directory containing protein sequences in fasta format for each genome will be created as well as a directory containing the corresponding nucleotide seqs.

Requirements: IPython, Prokka (http://www.vicbioinformatics.com/software.prokka.shtml)

Current Versions: IPython v 0.13.1, Prokka v 1.7

Usage: RunProkka.ipy [fasta directory]
