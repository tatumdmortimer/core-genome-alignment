#!/usr/bin/env python

import sys
import os
import argparse

# This script filters protein groups output from OrthoMCL and compares the core
# genomes of groups input by the user

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename


def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Compare core genomes')
    parser.add_argument("groups", help="OrthoMCL groups output", type=is_file)
    parser.add_argument("categories", help="File describing genome categories",
        type=is_file)
    parser.add_argument("-n", "--number", help="Number of genomes in a \
category allowed to not have a gene (default: 0)", type=int, default=0)
    parser.add_argument("-d", "--duplicates", help="Allow duplicate genes in \
one genome", action='store_true')
    return parser.parse_args()

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

def get_core_genomes(groupsDict, catDict, number, duplicates):
    """ Gets core genome for each category """ 
    coreDict = {}
    allDict = {} # holds genes present in category, but not shared in all
    for cat in catDict:
        coreDict[cat] = set()
        allDict[cat] = set()
    for group in groupsDict:
        genomeList = []
        proteinList = groupsDict[group]
        for protein in proteinList:
            ids = protein.split('|')
            genomeID = ids[0]
            genomeList.append(genomeID)
        genomeSet = set(genomeList)    # create set to check for duplicates
        for cat in catDict:
            if len(catDict[cat] - genomeSet) <= number:
                if duplicates:
                    coreDict[cat].add(group)
                elif len(genomeList) == len(genomeSet):
                    coreDict[cat].add(group)
            if len(catDict[cat] & genomeSet) > 0:
                allDict[cat].add(group)
    return coreDict, allDict

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

def get_unique_genes(coreDict, dupDict):
    uniqueDict = {}
    for group in coreDict:
        unique = coreDict[group].copy()
        for otherGroup in coreDict:
            if otherGroup != group:
                unique -= allDict[otherGroup]
        uniqueDict[group] = unique
    return uniqueDict

def write_unique_genomes(uniqueDict):
    for group in uniqueDict:
        outfile = open(group + "_uniqueGenes.txt", 'w')
        for gene in uniqueDict[group]:
            outfile.write(gene + '\n')
        outfile.close()

args = get_args()

groupsDict = read_groups_file(args.groups)
catDict = read_cat_file(args.categories)
coreDict, allDict = get_core_genomes(groupsDict, catDict, args.number,
args.duplicates)
write_core_genomes(groupsDict, catDict, coreDict)
uniqueDict = get_unique_genes(coreDict, allDict)
write_unique_genomes(uniqueDict)
