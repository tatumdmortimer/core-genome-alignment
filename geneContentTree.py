#!/usr/bin/env python

import sys
import os
import getopt
import re

sys.path.insert(1, "/opt/PepPrograms/biopython")
try:
    from Bio import Phylo
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
    from Bio.Phylo.TreeConstruction import _DistanceMatrix
except ImportError:
    print "oops, the import didn't work"
    sys.exit()

# This script uses the orthologous groups output by OrthoMCL to create a    
# distance matrix between strains. A neighbor joining tree is created from this
# distance matrix.

def get_arguments(argv):
    if len(argv) == 0:
        usage()
        sys.exit(2)
    groupsFile = None
    try:
        opts, args = getopt.getopt(argv, "g:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-g':
            groupsFile = arg
    return groupsFile

def usage():
    print "geneContentTree.py\n \
        -g <OrthoMCL groups file>"

def read_groups(groupsFile):
    print "Reading groups file"
    strainList = []
    strainDict = {}
    inFile = open(groupsFile, 'r')
    for line in inFile:
        line = line.strip()
        entries = line.split(':')
        groupName = entries[0]
        groupProteins = entries[1][1:].split(' ')
        for protein in groupProteins:
            strain = re.findall(r"[a-zA-Z0-9]+", protein)[1]
            if strain in strainDict:
                strainDict[strain].add(groupName)
            else:
                strainList.append(strain)
                strainDict[strain] = {groupName}
    inFile.close()
    print "Done reading groups file"
    return strainList, strainDict

def calc_dist(totalA, totalB, totalShared):
    distance = (totalA + totalB - 2.0*totalShared)/(totalA + totalB)
    print distance

def create_distance_matrix(strainList, strainDict):
    print "Calculating distance matrix"
    print strainList
    matrix = []
    for i in range(1, len(strainList) + 1):
        matrix.append([0]*i)
    dm = _DistanceMatrix(strainList, matrix)
    for a in range(len(strainList)):
        for b in range(a, len(strainList)):
            strA = strainList[a]
            strB = strainList[b]
            genA = strainDict[strA]
            genB = strainDict[strB]
            dm[strA, strB] = calc_dist(len(genA), len(genB), len(genA & genB))
    print "Done calculating distance matrix"
    return dm

def nj_tree(distanceMatrix):
    print "Constructing Neighbor Joining Tree"
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distanceMatrix)
    Phylo.write(tree, "geneContentTree.newick", "newick")
    print "Done constructing tree"

if None is get_arguments(sys.argv[1:]):
    usage()
    sys.exit(2)
else:
    groupsFile = get_arguments(sys.argv[1:])

strainList, strainDict = read_groups(groupsFile)
distanceMatrix = create_distance_matrix(strainList, strainDict)
nj_tree(distanceMatrix)
