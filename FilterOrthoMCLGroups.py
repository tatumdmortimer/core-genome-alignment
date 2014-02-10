#!/usr/bin/env python

import sys
import os
import getopt

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
        -n <number of genomes>"

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

def filter_groups(groupsDict, genomeNum):
    """ Looks through groups and checks if they have a single protein for 
    each genome. Returns a list of groups that meets this criteria"""
    keepList = []
    for group in groupsDict:
        genomeList = []
        proteinList = groupsDict[group]
        for protein in proteinList:
            ids = protein.split('|')
            genomeID = ids[0]
            genomeList.append(genomeID)
        genomeSet = set(genomeList)    # create set to check for duplicates
        if len(genomeList) == genomeNum and len(genomeSet) == genomeNum:
            keepList.append(group)
    return keepList

def write_file(inFileName, groupsDict, keepList):
    """ Writes groups in keep list with their proteins to a file."""
    outFileName = os.path.splitext(inFileName)[0] + '_filter.txt'
    outFile = open(outFileName, 'w')
    for group in keepList:
        outFile.write(group + ': ' + ' '.join(sorted(groupsDict[group])) + '\n')
    outFile.close()

if None in get_arguments(sys.argv[1:]):
    usage()
    sys.exit(2)
else:
    inFileName, genomeNum = get_arguments(sys.argv[1:])

groupsDict = read_groups_file(inFileName)
keepList = filter_groups(groupsDict, genomeNum)
write_file(inFileName, groupsDict, keepList)
