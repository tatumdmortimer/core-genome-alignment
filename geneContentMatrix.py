#!/usr/bin/env python

import sys
import os
import getopt
import re

# This script uses the orthologous groups output by OrthoMCL to create a    
# a matrix of gene content.

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
    groupList = []
    groupDict = {}
    inFile = open(groupsFile, 'r')
    for line in inFile:
        line = line.strip()
        entries = line.split(':')
        groupName = entries[0]
        groupList.append(groupName)
        groupDict[groupName] = []
        groupProteins = entries[1][1:].split(' ')
        for protein in groupProteins:
            strain = re.findall(r"[a-zA-Z0-9]+", protein)[0]
            groupDict[groupName].append(strain)
            if strain not in strainList:
                strainList.append(strain)
    inFile.close()
    print "Done reading groups file"
    return strainList, groupList, groupDict

def write_matrix(strainList, groupList, groupDict):
    print "Writing matrix file"
    outFile = open("geneContentMatrix.txt", "w")
    outFile.write("Group\t" + "\t".join(strainList) + "\n")
    for group in groupList:
        outFile.write(group)
        for strain in strainList:
            outFile.write("\t" + str(groupDict[group].count(strain)))
        outFile.write("\n")
    outFile.close()
     
if None is get_arguments(sys.argv[1:]):
    usage()
    sys.exit(2)
else:
    groupsFile = get_arguments(sys.argv[1:])

strainList, groupList, groupDict = read_groups(groupsFile)
write_matrix(strainList, groupList, groupDict)

