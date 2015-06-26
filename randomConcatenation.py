#!/usr/bin/python

import sys, os, argparse, random
from Bio import AlignIO
from Bio import SeqIO

################################################################################
# This script makes a core genome alignment from proteins in a groups file in
# the format output by OrthoMCL. FilterOrthoMCLGroups.py should be used first
# fasta files with sequences for proteins (AA or DNA) need to be in a directory 
# that is provided for the script.
################################################################################

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
            os.path.abspath(os.path.expanduser(values)))

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

def is_dir(dirname):
    """Checks if a path is a directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Randomly concatenate alignments')
    parser.add_argument("alignments", 
        help="Directory of alignments to choose from in fasta format",
        action=FullPaths, type=is_dir)
    parser.add_argument("-n", "--number", 
        help="Number of alignments to concatenate (default: 5)",
        type=int, default=5)
    parser.add_argument("-s", "--samples", 
        help="Number of concatenated alignments to create (default: 1000)",
        type=int, default=1000)
    return parser.parse_args()

def get_alignments(directory):
    alignments = listdir_fullpath(directory)
    firstAlign = AlignIO.read(alignments[0], "fasta")
    seqIDs = []
    for record in firstAlign:
        seqIDs.append(record.id.split('_')[0])
    return (alignments, seqIDs)

def concatenate_alignments(seqIDs, fileList, sample_num):
    """ Using Biopython, concatenate alignments """
    print "Concatenating gene alignments {0}".format(sample_num)
    inFileList = []
    for f in fileList:
        if os.path.isfile(f):
            inFileList.append(f)
        else:
            gene = f.split("_")[0]
            print ("%s does not exist" % gene)

    alignment = AlignIO.read(inFileList[0], 'fasta')
    num = len(alignment)
    for i in range(1, len(inFileList)):
        align = AlignIO.read(inFileList[i], 'fasta')
        if len(align) == num:
            alignment = alignment + align
        else:
            print inFileList[i].split("_")[0] + " doesn't have correct strain #"
    for i in range(len(seqIDs)):
        alignment[i].id = seqIDs[i]
    AlignIO.write(alignment, "alignment{0}.fasta".format(sample_num), 'fasta')

args = get_args()

alignments, seqIDs = get_alignments(args.alignments)
for i in range(args.samples):
    s = random.sample(alignments, args.number)
    concatenate_alignments(seqIDs, s, i)
