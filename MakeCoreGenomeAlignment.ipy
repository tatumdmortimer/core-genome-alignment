#!/usr/bin/env python

import sys
import os
import getopt

# This script makes a core genome alignment from proteins in a groups file in
# the format output by OrthoMCL. FilterOrthoMCLGroups.py should be used first
# fasta files with DNA sequences for proteins need to be in a directory that
# is provided for the script.

def get_arguments(argv):
    if len(argv) == 0:
        usage()
        sys.exit(2)
    inGroupsFile = None
    dnaSeqDirectory = None
    try:
        opts, args = getopt.getopt(argv, "g:d:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-d':
            dnaSeqDirectory = arg
        elif opt == '-g':
            inGroupsFile = arg
    return (inGroupsFile, numberOfGenomes)

def usage():
    print "MakeCoreGenomeAlignment.py\n \
        -g <input groups file>\n \
        -d <directory with fasta files of DNA sequences for proteins>"

def read_groups_file(inFileName):
    """ Read in groups file and create dictionary of group name and proteins in
    group"""
    inFile = open(inFileName, 'r')
    groupsDict = {}
    for line in inFile:
        line = line.strip()
        entries = line.split(':')
        groupName = entries[0]
        groupProteins = entries[1][1:].split(' ')
        groupsDict[groupName] = groupProteins
    inFile.close()
    return groupsDict

if None in get_arguments(sys.argv[1:]):
    usage()
    sys.exit(2)
else:
    inFileName, dnaDirectory = get_arguments(sys.argv[1:])

groupsDict = read_groups_file(inFileName)
