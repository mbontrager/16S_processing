#!/usr/bin/python

import os, sys, getopt, subprocess, glob, re
from parse_csv import csv_parse

############################################################
# Generate mothur commands and run mothur
#
# Input - path to demultiplexed .fasta files 
#
# Usage: Run from script directory:
# mothur_run.py -c <fasta_path/> -p <processors>
#
#
# Author: Martin Bontrager
############################################################


def main():
    proc = '1'
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hc:p:",["path=", "processors="])
    except getopt.GetoptError:
        print('mothur_run.py -c <fasta_path/> -p <processors>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('mothur_run.py -c <fasta_path/> -p <processors>')
            sys.exit()
        elif opt in ("-p", "--processors"):
            proc = str(arg)
        elif opt in ("-c", "--path"):
            path = arg
                   
    mothur(path, proc)
    
# Generate a make.contigs() mothur command with all samples
def gen_makegroups(path):
    
    os.chdir(path)
    full = glob.glob('*.fasta')
    files = []
    for f in full:
        if os.path.getsize(path + '/' + f) == 0:
            continue
        else:
            files.append(f)
    
    groupcmd = 'make.group(fasta='
    mergecmd = 'merge.files(input='
    for f in files:
        groupcmd = groupcmd + f + '-'
        mergecmd = mergecmd + f + '-'
    groupcmd = groupcmd.rstrip('-') + ', groups='
    mergecmd = mergecmd.rstrip('-') + ', output=allsamples.fasta)\n'
    for f in files:
        groupcmd = groupcmd + f.replace('.fasta', '') + '-'
    groupcmd = groupcmd.rstrip('-') + ')\n'
    mothurcmd = groupcmd + mergecmd
    return mothurcmd

# Generate initial mothur commands and run mothur batch file
def mothur(path, p):
    f = open('batch.mothur', 'r')
    path1 = path.replace('quality_filtered/', '')
    tmp = f.readlines()
    f.close()
    f = open((path1 + 'batch.mothur'), 'w')
    for i, line in enumerate(tmp):
        f.write(line)
        if i == 0:
            f.write('set.dir(input=' + path + ')\n')
            f.write(gen_makegroups(path))
            f.write('set.dir(input=' + path1 + 'mothur, output='+ path1 + 'mothur)\n')
            f.write('system(mv ' + path + 'mergegroups ' + path1 + 'mothur)\n')
            f.write('summary.seqs(fasta=allsamples.fasta, processors=' + p + ')\n')
    f.close()

    #os.chdir(path + 'mothur/')
#    run('mothur ../batch.mothur')

if __name__ == '__main__':
  main()
