#!/usr/bin/python

import sys, os, subprocess, shlex, glob, argparse
from datetime import datetime
from multiprocessing.dummy import Pool as ThreadPool

################################################################################
# This script runs the program OrthoMCL to put proteins from annotated genomes
# into orthologous groups. A sql database must be created and empty for this
# program to run correctly. Details of the database should be provided in the
# OrthoMCL config file. Additionally, the directory where the annotated proteins
# are located should be provided to this script. Do not use relative paths.
################################################################################

ORTHOMCL_PATH = "/opt/PepPrograms/orthomclSoftware-v2.0.9/bin/"

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
    parser = argparse.ArgumentParser(description='Run OrthoMCL pipeline')
    parser.add_argument("config", help="OrthoMCL config file", action=FullPaths,
        type=is_file)
    parser.add_argument("proteins",
        help="Directory with .faa files of proteins for each genome",
        action=FullPaths, type=is_dir)
    parser.add_argument("-t", "--threads",
        help="Number of threads to use (default: 1)",
        type=int, default=1)
    return parser.parse_args()

def call_with_log(cmd):
    """Calls a system command with the subprocess module. Redirects both stdout
    and stderr to a log file"""
    cmd = cmd.format(**(kvmap))

    logfile = open(orthomcl_wd + current_datetime+".log", "a+")
    logfile.write("Executing command: " + cmd + "\n")
    logfile.flush()
    ret = subprocess.call(shlex.split(cmd), stdout=logfile, stderr=logfile)
    if(ret != 0):
        print("Pipeline did not complete successfully. \n Command : \n\n" +
            cmd + "\n\n returned with non-zero code: " + str(ret))
        logfile.write("Pipeline did not complete successfully.\nCommand :\n\n" +
            cmd + "\n\n returned with non-zero code: " + str(ret))
        logfile.close()
        sys.exit(-1)
    logfile.close()

def install_schema(config):
    print("Installing schema")
    call_with_log(ORTHOMCL_PATH + "orthomclInstallSchema %s log.sql.txt" %
        config)

def compliant_fasta(directory):
    def edit_fasta(fasta_file):
        prefix = os.path.splitext(os.path.basename(fasta_file))[0]
        if len(prefix) > 4:
            prefix = prefix[-4:]
        call_with_log(ORTHOMCL_PATH +
        "orthomclAdjustFasta %s %s 1" % (prefix, fasta_file))
        return prefix

    print "Making fasta files compliant"
    pool = ThreadPool(args.threads)
    fasta_files = listdir_fullpath(directory)
    prefixes = pool.map(edit_fasta, fasta_files)
    pool.close()
    pool.join()
    dups = set([x for x in prefixes if prefixes.count(x) > 1])
    if len(dups) > 0:
        print ("WARNING: duplicate prefixes")

def filter_fasta():
    print "Filtering fasta files"
    call_with_log(ORTHOMCL_PATH + "orthomclFilterFasta compliantFASTA/ 10 20")

def call_blast_parser(cmd, filename):
    cmd = cmd.format(**(kvmap))
    outfile = open(filename, "w+")
    logfile = open(orthomcl_wd + current_datetime + ".log", "a+")
    logfile.write("Executing command: " + cmd + "\n")
    logfile.flush()
    ret = subprocess.call(shlex.split(cmd), stdout=outfile, stderr=logfile)
    logfile.close()
    outfile.close()

def allVSall_blast(config):
    print "Performing all vs. all blast"
    call_with_log("makeblastdb -in ../goodProteins.fasta -dbtype\
    prot -out coreDB")
    call_with_log("blastp -db coreDB -query ../goodProteins.fasta \
-outfmt 6 -out allVSall.tsv -num_threads %i" % args.threads)
    call_blast_parser(ORTHOMCL_PATH +
"orthomclBlastParser allVSall.tsv ../compliantFASTA/", "similarSequences.txt")
    call_with_log(ORTHOMCL_PATH + "orthomclLoadBlast %s similarSequences.txt"
         % config)

def pairs(config):
    print "Performing pairs step"
    call_with_log(ORTHOMCL_PATH + "orthomclPairs %s pairs.log cleanup=yes" %
        config)
    call_with_log(ORTHOMCL_PATH + "orthomclDumpPairsFiles %s" % config)

def call_groups_log(cmd, infilename, outfilename):
    cmd = cmd.format(**(kvmap))
    outfile = open(outfilename, "w+")
    infile = open(infilename, "r")
    logfile = open(orthomcl_wd + current_datetime + ".log", "a+")
    logfile.write("Executing command: " + cmd + "\n")
    logfile.flush()
    ret = subprocess.call(shlex.split(cmd), stdin=infile,
        stdout=outfile, stderr=logfile)
    logfile.close()
    outfile.close()


def mcl():
    print "Running mcl"
    call_with_log("mcl mclInput --abc -I 1.5 -o mclOutput")
    call_groups_log(ORTHOMCL_PATH + "orthomclMclToGroups core 1000",
        "mclOutput", "groups.txt")


args = get_args()
current_datetime = datetime.today().strftime("%d-%m-%Y-%H%M")
orthomcl_wd = os.getcwd() + "/"
kvmap = {'projectname':'orthomcl'}

install_schema(args.config)
try:
    os.makedirs("compliantFASTA")
except OSError:
    if not os.path.isdir("compliantFASTA"):
        raise
os.chdir('compliantFASTA')
compliant_fasta(args.proteins)
os.chdir('..')
filter_fasta()
try:
    os.makedirs("blast")
except OSError:
    if not os.path.isdir("blast"):
        raise
os.chdir('blast')
allVSall_blast(args.config)
os.chdir('..')
pairs(args.config)
mcl()

