#!/usr/bin/python

import os, sys, getopt, subprocess, glob
from parse_csv import csv_parse

############################################################
# Demultiplex 338F/806R 16S MiSeq Reads with internal barcode
#
# Input - Paired .fastq reads w/ internal barcodes and a 
# barcode mapping file in .csv format. See README for details
#
# Usage: demux_samples.py -f forward.fastq -r reverse.fastq -c map.csv
#
# Author: Martin Bontrager
############################################################

def main():
    f338 = ''
    r806 = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hf:r:c:",["forward=","reverse=", "csvfile="])
    except getopt.GetoptError:
        print('demux_samples.py -1 <forward_read.fastq> -2 <reverse_read.fastq> -c <map.csv>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('demux_samples.py -1 <forward_read.fastq> -2 <reverse_read.fastq> -c <map.csv>')
            sys.exit()
        elif opt in ("-f", "--forward"):
            f338 = arg
        elif opt in ("-r", "--reverse"):
            r806 = arg
        elif opt in ("-c", "--csvfile"):
            csvfile = arg
    
    path = os.path.dirname(f338) + '/'
    csv_parse(csvfile, (path + 'barcodes.fil'), (path + 'samples.txt'))
    subprocess.call(['mkdir', (path + 'demux')])
    subprocess.call(['mkdir', (path + 'reports')])
    subprocess.call(['mkdir', (path + 'reports/FLASH_histograms')])
    subprocess.call(['mkdir', (path + 'overlapped')])
    demux(f338, r806, path)
    dpath = path + 'demux/'

    for i in get_samples(path + 'samples.txt'):
        trim_barcodes((dpath + i + '_338F.fastq'), (dpath + i + '_806R.fastq'))
        overlap(i, (dpath + i + '_338F_bctrimmed.fastq'), (dpath + i + '_806R_bctrimmed.fastq'))
    p = subprocess.Popen('mv ' + (path + 'overlapped/*.hist ') + 
                         (path + 'reports/FLASH_histograms/'), shell=True)
    p.communicate()
    subprocess.Popen('rm ' + (path + 'overlapped/*.histogram'), shell=True)
    p.communicate()

# Simplify running bash commands
def run(cmd):
    p = subprocess.Popen(cmd, shell=True)
    os.waitpid(p.pid, 0)

# Demultiplex the samples
def demux(forward, reverse, path):
    cmd = ('../tools/fastq-multx -x -b -B ' + path + 'barcodes.fil ' + forward + ' ' + reverse + ' -o ' + 
           path + 'demux/%_338F.fastq ' + path + 'demux/%_806R.fastq -m 1 > ' + path + 
           'reports/demux_log.txt')
    run(cmd)

# Trim 12bp barcodes from each sequence read
def trim_barcodes(forward, reverse):
    newfor = forward.replace('.fastq', '') + '_bctrimmed.fastq'
    newrev = reverse.replace('.fastq', '') + '_bctrimmed.fastq'
    cmd = ('../tools/seqtk/seqtk trimfq -b 12 ' + forward + ' > ' + newfor)
    run(cmd)
    subprocess.call(['rm', forward])
    cmd = ('../tools/seqtk/seqtk trimfq -b 12 ' + reverse + ' > ' + newrev)
    run(cmd)
    subprocess.call(['rm', reverse])

# Get a list of samples from the 'samples.txt' file generated via csv_parse
def get_samples(samples):
    s = []
    with open(samples) as f:
        for l in f:
            s.append(l.strip())
    return s

# Overlap reads using flash
def overlap(sample, forward, reverse):
    d = os.path.dirname(forward).replace('/demux', '/overlapped')
    cmd = ('../tools/FLASH-1.2.8/flash ' + forward + ' ' +  reverse + 
           ' -r 288 -f 429 -s 18 -d '+ d + ' -o ' + sample + ' > ' + 
           d.replace('overlapped', 'reports/') + sample + '.flash.log')
    run(cmd)    
    

if __name__ == '__main__':
  main()
