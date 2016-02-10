#!/usr/bin/python

import csv, os, sys, getopt, subprocess, glob

############################################################
# Parse IGS .csv file with barcode names, sample ID, etc. 
# Generate a "barcodes.fil" for fastq-multx demux
#
# Specific to output from the Ravel lab MiSeq 16S pipe
#
# Usage: python parse_csv.py -i infile.csv -o outfile
#
# Author: Martin Bontrager
############################################################

def main():
    csv_file = ''
    outfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["input=", "output="])
    except getopt.GetoptError:
        print('parse_csv.py -i input.csv -o outfile')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('parse_csv.py -i input.csv -o outfile')
            sys.exit()
        elif opt in ("-i", "--input"):
            csv_file = arg
        elif opt in ("-o", "--output"):
            outfile = arg
    csv_parse(csv_file, outfile)

# Parse the .csv barcodes file
def csv_parse(file, outfile, samples):
    output = open(outfile, 'w')
    samples = open(samples, 'w')
    output.write("# sample\tforward\treverse\tstyle\n")
    with open(file) as f:
        reader = csv.DictReader(f)
        for r in reader:
            outstring = (str(r['sample_name']) + '.' + str(r['row']) + str(r['well']) + '\t' 
                         + str(r['F_barcode']) + '-' + str(r['R_Barcode']) + '\tTruSeq DNA\n')
            output.write(outstring)
            samples.write(str(r['sample_name']) + '.' + str(r['row']) + str(r['well']) + '\n')
    output.close()
    samples.close()
if __name__ == "__main__":
    main()
