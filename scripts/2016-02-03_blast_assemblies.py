#!/usr/bin/python

import os, sys, getopt, glob, subprocess, shutil
from Bio import SeqIO

############################################################
# Align de novo assemblies against the reference copepod genome
# and process the output into a format for R analysis
#
# Usage: Run from script directory:
# python3 -d blast_database -p path_to_assemblies
#
# Author: Martin Bontrager
############################################################

def main():
    path = ''
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],'hp:d:',['path=', 'database='])
    except getopt.GetoptError:
        print('Fix your input')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('-d blast_db -p path_to_assemblies')
            sys.exit()
        elif opt in ("-p", "--path"):
            path = str(arg)
        elif opt in ("-d", "--database"):
            db = arg
            
    os.chdir(path)
    sample_list = glob.glob('*.fna')
    print(sample_list)
    print(path)
    print(db)
    

# Filter assemblies to only those contigs greater than n base pairs
def filter_by_length(input_file, length):
    infile = open(input_file, 'rU')
    ofname = (input_file.replace('__FINAL_ASSEMBLY.consolidatedContigs.fna',
                             '_3000.fna'))
    out = open(ofname, 'w')

    for i in SeqIO.parse(infile, 'fasta'):
	sequence = i.seq
	if len(sequence) > length:
		SeqIO.write(i, out, 'fasta')

    infile.close()
    out.close()
        
# Simplify running bash commands
def run(cmd):
    p = subprocess.Popen(cmd, shell=True)
    p.communicate()
    
if __name__ == "__main__":
    main()

    
