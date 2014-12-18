#!/usr/bin/python

import sys, os, argparse, subprocess, shlex, glob
from datetime import datetime
from multiprocessing.dummy import Pool as ThreadPool
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
    parser = argparse.ArgumentParser(description='Align and concatenate core \
genes from OrthoMCL groups file')
    parser.add_argument("groups", help="OrthoMCL groups file", action=FullPaths,
        type=is_file)
    parser.add_argument("sequences", 
        help="Directory with .fasta files of proteins (AA or DNA) for genomes",
        action=FullPaths, type=is_dir)
    parser.add_argument("-t", "--threads", 
        help="Number of threads to use (default: 1)",
        type=int, default=1)
    parser.add_argument("-c", "--codon",
        help="Align nucleotides using translation into amino acids",
        action='store_true')
    return parser.parse_args()

def call_with_log_mafft(cmd, outfilename):
    """Calls a system command with the subprocess module. Redirects both stdout
    and stderr to a log file"""
    cmd = cmd.format(**(kvmap))
    outfile = open(outfilename, "w+")
    logfile = open(align_wd + current_datetime+".log", "a+")
    logfile.write("Executing command: " + cmd + "\n")
    logfile.flush()
    ret = subprocess.call(shlex.split(cmd), stdout=outfile, stderr=logfile)
    if(ret != 0):
        print("Pipeline did not complete successfully. \n Command : \n\n" + 
            cmd + "\n\n returned with non-zero code: " + str(ret))
    logfile.close()
    outfile.close()

def call_with_log(cmd):
    """Calls a system command with the subprocess module. Redirects both stdout
    and stderr to a log file"""
    cmd = cmd.format(**(kvmap))
    logfile = open(align_wd + current_datetime+".log", "a+")
    logfile.write("Executing command: " + cmd + "\n")
    logfile.flush()
    ret = subprocess.call(shlex.split(cmd), stdout=logfile, stderr=logfile)
    if(ret != 0):
        print("Pipeline did not complete successfully. \n Command : \n\n" + 
            cmd + "\n\n returned with non-zero code: " + str(ret))
    logfile.close()

def read_groups_file(inFileName):
    """ Read in groups file and create dictionary of group name and proteins in
    group"""
    print "Reading groups file"
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

def make_unaligned_fasta(dnaDirectory, groupsDict):
    """ Reads through files in provided directory to find gene sequences that
    match the proteins in the groups dictionary"""
    print "Collecting core genes"
    def make_fasta(group):
        proteins = groupsDict[group]
        out = open('proteinAlignments/' + group + '.fasta', 'w')
        records = []
        seqIDs = []
        for protein in proteins:
            seqID = protein.split('|')[0]
            seqIDs.append(seqID)
            protein = protein.split('|')[1]
            records.append(seqRecordDict[protein])
        SeqIO.write(records, out, 'fasta')
        return seqIDs

    try:
        os.makedirs("proteinAlignments")
    except OSError:
        if not os.path.isdir("proteinAlignments"):
            raise
    files = listdir_fullpath(dnaDirectory)
    seqRecordDict = {}
    seqIDs = []
    for f in files:
        handle = open(f, 'r')
        for record in SeqIO.parse(handle, 'fasta'):
            seqRecordDict[record.id] = record
    pool = ThreadPool(args.threads)
    seqIDs = pool.map(make_fasta, groupsDict.keys())
    pool.close()
    pool.join()
    return seqIDs[0]


def align_gene_sequences(groupFastaDirectory):
    """ Use MAFFT to align gene sequences"""
    print "Aligning core genes"
    def run_mafft(infile):
        outfile = os.path.splitext(infile)[0] + '_align.fasta'
        call_with_log_mafft("mafft --globalpair --maxiterate 1000 --thread 1 %s"
            % infile, outfile)
        return outfile
    def run_translatorX(infile):
        outfile = os.path.splitext(infile)[0] 
        call_with_log("/opt/PepPrograms/translatorx_vLocal.pl -i %s -o %s -p F"
            % (infile, outfile))
        return outfile + ".nt_ali.fasta"
    files = listdir_fullpath(groupFastaDirectory)
    pool = ThreadPool(args.threads)
    if args.codon:
        outfileList = pool.map(run_translatorX, files)
    else:
        outfileList = pool.map(run_mafft, files)
    pool.close()
    pool.join()
    return outfileList

def trim_alignments(fileList):
    """ Use trimAl to trim gaps from input alignments"""
    print "Trimming alignments"
    def run_trimal(infile):
        outfile = os.path.splitext(infile)[0] + '_trim.fasta'
        call_with_log("/opt/PepPrograms/trimal/source/trimal -in %s -out %s \
-automated1" % (infile, outfile))
        return outfile
    pool = ThreadPool(args.threads) 
    outfileList = pool.map(run_trimal, fileList)
    pool.close()
    pool.join()
    return outfileList

def concatenate_alignments(seqIDs, fileList):
    """ Using Biopython, concatenate alignments produced by MAFFT"""
    print "Concatenating gene alignments"
    inFileList = []
    for f in fileList:
        if os.path.isfile(f):
            inFileList.append(f)
        else:
            gene = f.split("_")[0]
            print ("%s is not included in core alignment" % gene)

    coreAlignment = AlignIO.read(inFileList[0], 'fasta')
    num = len(coreAlignment)
    for i in range(1, len(inFileList)):
        align = AlignIO.read(inFileList[i], 'fasta')
        if len(align) == num:
            coreAlignment = coreAlignment + align
        else:
            print inFileList[i].split("_")[0] + " doesn't have correct strain #"
    for i in range(len(seqIDs)):
        coreAlignment[i].id = seqIDs[i]
    AlignIO.write(coreAlignment, "core_alignment.fasta", 'fasta')

args = get_args()
current_datetime = datetime.today().strftime("%d-%m-%Y-%H%M")
align_wd = os.getcwd() + "/"
kvmap = {'projectname':'coreAlignment'}

groupsDict = read_groups_file(args.groups)
seqIDs = make_unaligned_fasta(args.sequences, groupsDict)
alignList = align_gene_sequences('proteinAlignments/')
if not args.codon:
    trimList = trim_alignments(alignList)
    concatenate_alignments(seqIDs, trimList)
else:
    concatenate_alignments(seqIDs, alignList)
