#!/usr/bin/env python

import sys
import os
import getopt
import glob

# This script reads in a fasta alignment or a directory of alignments.
# Alignments in directory must end in .fasta.
# Alignments must also be in frame coding sequence. 

def get_arguments(argv):
    if len(argv) == 0:
        usage()
        sys.exit(2)
    alignmentFile = None
    alignmentDirectory = None
    try:
        opts, args = getopt.getopt(argv, "a:d:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-a':
            alignmentFile = arg
        elif opt == '-d':
            alignmentDirectory = arg
    return (alignmentFile, alignmentDirectory)

def usage():
    print "editStrainNames.py\n \
        -a <fasta alignment>\n \
        -d <directory of fasta alignments>"

def editStrains(align):
    infile = open(align, 'r')
    outfile = open(os.path.splitext(align)[0] + '_strainsEdited.fasta', 'w')
    for line in infile:
        line = line.strip()
        if line[0] == '>':
            line = line.split('_')
            outfile.write(line[0] + '\n')
        else:
            outfile.write(line + '\n')
    infile.close()
    outfile.close()
       
    
alignment, directory = get_arguments(sys.argv[1:])
alignDict = {}
# Check if alignment or directory was given and calculate stats accordingly
if alignment is None:
    if directory is None:
        usage()
        sys.exit()
    else:
        for i,align in enumerate(glob.glob(directory + '*.fasta')):
            editStrains(align)

elif alignment is not None:
    if directory is not None:
        print "Must only input an alignment or a directory"
        usage()
        sys.exit()
    else:
        editStrains(alignment)

