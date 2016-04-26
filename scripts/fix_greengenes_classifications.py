#!/usr/bin/python

import os, sys, getopt, re

############################################################
# Rename greengenes classifications
#
# Input - mothur cons.taxonomy file
#
# Usage: Run from script directory:
# fix_greengenes_classifications.py -i mothur.cons.taxonomy
#
# Change as needed for each run of samples
#
# Author: Martin Bontrager
############################################################



def main():

    repl = (('k__', ''), ('p__', ''), ('c__', ''), ('o__', ''), ('f__', ''), ('g__', ''), ('s__', ''))
    infile = ''

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:",["input="])
    except getopt.GetoptError:
        print('fix_greengenes_classifications.py -i mothur.cons.taxonomy')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('fix_greengenes_classifications.py -i mothur.cons.taxonomy')
            sys.exit()
        elif opt in ("-i", "--input"):
            infile = arg
    
    change_string(infile, repl)

# Parse the .cons.taxonomy file to replace 'level__' tags and extraneous classification
def change_string(infile, repl):
    output = open((infile.replace('.cons.taxonomy', '') + '_fixed.cons.taxonomy'), 'wb')
    
    with open(infile, 'rb') as f:
        for r in f:
            s = re.sub(r'(\tk__.*;.*;.*;.*;.*;.*;.*;).*;', r'\1', r)
            v = reduce(lambda a, kv: a.replace(*kv), repl, s)
            output.write(v)
    output.close()

if __name__ == "__main__":
    main()
