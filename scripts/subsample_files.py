#!/usr/bin/python

import os, sys, getopt, glob, subprocess, shutil

############################################################
# Use mothur to subsample large fasta files
#
# Author: Martin Bontrager
############################################################

def main():
    path = ''
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],'hp:',['path='])
    except getopt.GetoptError:
        print('subsample_files -p path/to/fasta/dir')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('subsample_files -p path/to/fasta/dir')
            sys.exit()
        elif opt in ('-p', '--path'):
            path = arg
    
    os.chdir(path)
    sample_list = glob.glob('*.fasta')
    
    for f in sample_list:
        cmd = ('mothur135 "#sub.sample(fasta=' + f + ')"')
        run(cmd)
    
        
# Simplify running bash commands
def run(cmd):
    p = subprocess.Popen(cmd, shell=True)
    p.communicate()
    
if __name__ == "__main__":
    main()
