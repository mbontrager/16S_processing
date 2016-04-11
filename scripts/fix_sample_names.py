#!/usr/bin/python

import os, sys, getopt

############################################################
# Rename samples from IGS csvfile output
#
# Input - CSV file 
#
# Usage: Run from script directory:
# fix_sample_names.py -c map.csv
#
# Change as needed for each run of samples
#
# Author: Martin Bontrager
############################################################

def main():

    repl = (('PONTE-1', 'POE1'), ('PONTE-2', 'POE2'), ('MILE_1U', 'MIE'),
            ('_1', ''), ('_2', ''), ('_3', ''), ('_4', ''),
            ('_5', ''), ('_6', ''), ('-2', ''), ('-1', ''), ('LOUW01', 'LOW'),
            ('SJPE', 'SJE'), ('BRADE', 'BRE'), ('V1W01', 'V1W'),
            ('LOUE', 'LOE'), ('BRADW', 'BRW'), ('MILW1', 'MIW'), ('VIE', 'V1E'),
            ('W-', 'W'), ('AE-', 'AE'), ('W(-)', 'W'), ('AE(-)', 'AE'))
    csv_file = ''

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hc:",["csvfile="])
    except getopt.GetoptError:
        print('fix_sample_names.py -c map.csv')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('fix_sample_names.py -c map.csv')
            sys.exit()
        elif opt in ("-c", "--csvfile"):
            csv_file = arg
    
    change_names(csv_file, repl)

# Parse the .csv barcodes file
def change_names(csvfile, repl):
    output = open((csvfile.replace('.csv', '') + '_fixed.csv'), 'wb')
    
    with open(csvfile, 'rb') as f:
        for r in f:
            v = reduce(lambda a, kv: a.replace(*kv), repl, r)
            output.write(v)
    output.close()

if __name__ == "__main__":
    main()

