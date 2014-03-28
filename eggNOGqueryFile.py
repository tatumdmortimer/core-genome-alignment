#!/usr/bin/env python

import sys
import os
import getopt
import glob
from Bio import SeqIO

# This script takes the first protein from each file in a directory of unaligned
# fasta files and changes the protein ID to the name of the file. All protein  
# sequences are written to an output file that should be used as the query for
# eggNOGblast.py.                                                   

def usage():
    print "FilterOrthoMCLGroups.py <directory of .fasta protein sequences>"

if len(sys.argv) != 2:
    usage()
    sys.exit()

protDir = sys.argv[1]
outfile = open("groupSeq.fasta", 'w')
proteins = []

for fasta in glob.glob(protDir + '*.fasta'):
    first_record = next(SeqIO.parse(fasta, "fasta"))
    first_record.id = os.path.splitext(fasta)[0].split('/')[-1]
    proteins.append(first_record)

SeqIO.write(proteins, outfile, 'fasta')
outfile.close()
