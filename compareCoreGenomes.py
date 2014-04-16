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
    genomeCatFile = None
    try:
        opts, args = getopt.getopt(argv, "g:c:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-c':
            genomeCatFile = arg
        elif opt == '-g':
            inGroupsFile = arg
    return (inGroupsFile, genomeCatFile)

def usage():
    print "compareCoreGenomes.py\n \
        -g <input groups file>\n \
        -c <input genome categories file"

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

def read_cat_file(genomeCatFile):
    """ Read in genome categories and create dictionary of category name and 
    genomes in that category"""
    inFile = open(genomeCatFile, 'r')
    catDict = {}
    for line in inFile:
        line = line.strip()
        entries = line.split()
        genome = entries[0]
        cat = entries[1]
        if cat in catDict:
            catDict[cat].add(genome)
        else:
            catDict[cat] = {genome}
    inFile.close()
    return catDict

def get_core_genomes(groupsDict, catDict):
    """ Gets core genome for each category """ 
    coreDict = {}
    for cat in catDict:
        coreDict[cat] = set()
    for group in groupsDict:
        genomeList = []
        proteinList = groupsDict[group]
        for protein in proteinList:
            ids = protein.split('|')
            genomeID = ids[0]
            genomeList.append(genomeID)
        genomeSet = set(genomeList)    # create set to check for duplicates
        if len(genomeList) == len(genomeSet):
            for cat in catDict:
                if catDict[cat].issubset(genomeSet):
                    coreDict[cat].add(group)
    return coreDict

def write_core_genomes(groupsDict, catDict, coreDict):
    for group in coreDict:
        outfile = open(group + "_coreGenes.txt", 'w')
        for gene in coreDict[group]:
            proteinList = groupsDict[gene]
            proteinsInCat = []
            for protein in proteinList:
                ids = protein.split('|')
                if ids[0] in catDict[group]:
                    proteinsInCat.append(protein)
            outfile.write(gene + ': ' + ' '.join(sorted(proteinsInCat)) + '\n')
        outfile.close()

def get_unique_genes(coreDict):
    uniqueDict = {}
    for group in coreDict:
        unique = coreDict[group]
        for otherGroup in coreDict:
            if otherGroup != group:
                unique -= coreDict[otherGroup]
        uniqueDict[group] = unique
    return uniqueDict

def write_unique_genomes(uniqueDict):
    for group in uniqueDict:
        outfile = open(group + "_uniqueGenes.txt", 'w')
        for gene in uniqueDict[group]:
            outfile.write(gene + '\n')
        outfile.close()

if None in get_arguments(sys.argv[1:]):
    usage()
    sys.exit(2)
else:
    inFileName, genomeCatFile = get_arguments(sys.argv[1:])

groupsDict = read_groups_file(inFileName)
catDict = read_cat_file(genomeCatFile)
coreDict = get_core_genomes(groupsDict, catDict)
write_core_genomes(groupsDict, catDict, coreDict)
uniqueDict = get_unique_genes(coreDict)
write_unique_genomes(uniqueDict)
