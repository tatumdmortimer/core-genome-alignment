#!/usr/bin/env python

import argparse
from Bio import Entrez
from Bio import SeqIO


# This script downloads genomes from NCBI using accession number

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Download genome')
    parser.add_argument("accession", help="NCBI Accession Number")
    parser.add_argument("email", help="E-mail address")
    return parser.parse_args()

args = get_args()

Entrez.email = args.email
handle = Entrez.efetch(db="nucleotide", id=args.accession, 
    rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
handle.close()

name = record.description.split(",")[0].split(" ")[1:]
genus = name[0]
species = name[1]
strain = name[2]
print genus, species, strain

SeqIO.write(record, "_".join(name) + ".fasta", "fasta")
