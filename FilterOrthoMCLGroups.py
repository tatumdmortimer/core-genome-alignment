#!/usr/bin/env python

import sys
import re
import getopt
from operator import itemgetter
import random

# This script filters protein groups output from OrthoMCL and returns those that
# contain only one protein per genome and have all genomes represented

def get_arguments(argv):
    if len(argv) == 0:
        usage()
        sys.exit(2)
    inGroupsFile = None
    numberOfGenomes = None
    try:
        opts, args = getopt.getopt(argv, "g:n:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-n':
            numberOfGenomes = int(arg)
        elif opt == '-g':
            inGroupsFile = arg
    return (inGroupsFile, numberOfGenomes)

def usage():
    print "FilterOrthoMCLGroups.py\n \
        -g <input groups file>\n \
        -o <number of genomes>"

if None in get_arguments(sys.argv[1:]):
    usage()
    sys.exit(2)
else:
    inFile, genomeNum = get_arguments(sys.argv[1:])
