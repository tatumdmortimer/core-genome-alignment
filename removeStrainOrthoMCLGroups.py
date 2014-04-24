#!/usr/bin/env python

import sys
import os
import getopt

# This script filters protein groups output from OrthoMCL and compares the core
# genomes of groups input by the user

def get_arguments(argv):
    if len(argv) == 0:
        usage()
        sys.exit(2)
    inGroupsFile = None
    remove = None
    try:
        opts, args = getopt.getopt(argv, "g:r:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-g':
            inGroupsFile = arg
        elif opt == '-r':
            remove = arg
    return inGroupsFile, remove

def usage():
    print "compareCoreGenomes.py\n \
        -g <input groups file>\n \
        -r <strain to remove from groups>"

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

def write_groups(groupsDict, remove):
    outfile = open("groups_no" + remove + '.txt', 'w')
    for group in sorted(groupsDict.keys()):
        proteinList = groupsDict[group]
        goodProteinList = []
        for protein in proteinList:
            if protein.split('|')[0] != remove:
                 goodProteinList.append(protein)
        outfile.write(group + ': ' + ' '.join(sorted(goodProteinList)) + '\n')
    outfile.close()

if None in get_arguments(sys.argv[1:]):
    usage()
    sys.exit(2)
else:
    inFileName, remove = get_arguments(sys.argv[1:])

groupsDict = read_groups_file(inFileName)
write_groups(groupsDict, remove)
