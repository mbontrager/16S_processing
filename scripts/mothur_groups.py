#!/usr/bin/python

import os, sys, getopt, subprocess, glob, re
from parse_csv import csv_parse

############################################################
# Simply generate mothur group commands for further edits
#
# Input - file with group names. One per line.
#
# Usage: Run from script directory:
# mothur_groups.py -f <file.txt>
#
#
# Author: Martin Bontrager
############################################################


def main():
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hf:",["file="])
    except getopt.GetoptError:
        print('mothur_groups.py -f <file.txt>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('mothur_groups.py -f <file.txt>')
            sys.exit()
        elif opt in ("-f", "--file"):
            f = arg
                   
    gen_groups(f)
    
# Generate a make.contigs() mothur command with all samples
def gen_groups(file):
    
    a = os.path.dirname(os.path.realpath(file))
    f = open((a + '/batch.mothur'), 'w')
    
    groupcmd = 'mothurcmd(groups='

    with open(file) as b: 
        for line in b:
            groupcmd = groupcmd + line.strip() + '-'
    
    groupcmd = groupcmd.rstrip('-') + ')\n'
    
    f.write(groupcmd)

if __name__ == '__main__':
  main()
