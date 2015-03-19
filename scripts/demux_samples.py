#!/usr/bin/python

import os, sys, getopt, subprocess, glob
from parse_csv import csv_parse

############################################################
# Demultiplex 338F/806R 16S MiSeq Reads with internal barcode
#
# Input - Paired .fastq reads w/ internal barcodes and a 
# barcode mapping file in .csv format. See README for details
#
# Usage: Run from script directory:
# demux_samples.py -f forward.fastq -r reverse.fastq -c map.csv
#
# Use the '-k' flag to keep intermediate seq files, otherwise 
# they will be deleted: log files will always be kept.
#
# Author: Martin Bontrager
############################################################

def main():
    forw = ''
    rev = ''
    keep = False

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hkf:r:c:",["forward=","reverse=", "csvfile="])
    except getopt.GetoptError:
        print('demux_samples.py -1 <forward_read.fastq> -2 <reverse_read.fastq> -c <map.csv>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('demux_samples.py -1 <forward_read.fastq> -2 <reverse_read.fastq> -c <map.csv>')
            sys.exit()
        if opt == '-k':
            keep = True
        elif opt in ("-f", "--forward"):
            forw = arg
        elif opt in ("-r", "--reverse"):
            rev = arg
        elif opt in ("-c", "--csvfile"):
            csvfile = arg
    
    path = os.path.dirname(forw) + '/'
    dpath = path + 'demux/'
    subprocess.call(['mkdir', (path + 'demux')])
    subprocess.call(['mkdir', (path + 'logs')])
    subprocess.call(['mkdir', (path + 'logs/FLASH_histograms')])
    subprocess.call(['mkdir', (path + 'logs/overlap')])
    subprocess.call(['mkdir', (path + 'logs/primer_trimming')])
    subprocess.call(['mkdir', (path + 'overlapped')])
    subprocess.call(['mkdir', (path + 'trimmed')])
    subprocess.call(['mkdir', (path + 'quality_filtered')])
    
    change_names(csvfile)
    csvfile = csvfile.replace('.csv', '_fixed.csv')
    csv_parse(csvfile, (path + 'barcodes.fil'), (path + 'samples.txt'))
    demux(forw, rev, path)
    
    for i in get_samples(path + 'samples.txt'):
        if check_files((path + 'demux/'), i, '_338F.fastq') == 0 :
           continue        
        trim_barcodes((dpath + i + '_338F.fastq'), (dpath + i + '_806R.fastq'))
        overlap_usearch(i, (dpath + i + '_338F_bctrimmed.fastq'), (dpath + i + '_806R_bctrimmed.fastq'))
        if check_files((path + 'overlapped/'), i, '.extendedFrags.fastq') == 0 :
           continue     
        trim_primers((path + 'overlapped/'), i)
        qc((path + 'trimmed/'), i)
    
    gen_makecontigs((path + 'quality_filtered'))
    # p = subprocess.Popen('mv ' + (path + 'overlapped/*.hist ') + 
    #                      (path + 'logs/FLASH_histograms/'), shell=True)
    # p.communicate()
    # subprocess.Popen('rm ' + (path + 'overlapped/*.histogram'), shell=True)
    # p.communicate()

# Simplify running bash commands
def run(cmd):
    p = subprocess.Popen(cmd, shell=True)
    p.communicate()

# Change poor file names with external script (modify as needed)
def change_names(csvfile):
    cmd = ('python fix_sample_names.py -c ' + csvfile)
    run(cmd)

# Demultiplex the samples
def demux(forward, reverse, path):
    cmd = ('../tools/fastq-multx -x -b -B ' + path + 'barcodes.fil ' + forward + ' ' + reverse + ' -o ' + 
           path + 'demux/%_338F.fastq ' + path + 'demux/%_806R.fastq -m 1 > ' + path + 
           'logs/demux_log.txt')
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
def overlap_flash(sample, forward, reverse):
    d = os.path.dirname(forward).replace('/demux', '/overlapped')
    cmd = ('../tools/FLASH-1.2.8/flash ' + forward + ' ' +  reverse + 
           ' -r 288 -f 429 -s 18 -d '+ d + ' -o ' + sample + ' > ' + 
           d.replace('overlapped', 'logs/') + sample + '.overlap.log')
    run(cmd)

# Overlap reads using usearch8 (preferred)
def overlap_usearch(sample, forward, reverse):
    d = os.path.dirname(forward).replace('/demux', '/overlapped/')
    cmd = ('../tools/usearch8 -fastq_mergepairs ' + forward + ' -reverse ' + 
           reverse + ' -fastq_minovlen 100 -fastqout '+ d + sample + 
           '.extendedFrags.fastq -log '+ d.replace('/overlapped/', '/logs/overlap/') +
           sample + '.overlap.log -quiet')
    run(cmd)    

# Trim 338F/806R primers with tagcleaner
def trim_primers(path, sample):
    cmd = ('../tools/tagcleaner-standalone-0.16/tagcleaner.pl -fastq ' + path + 
           sample + '.extendedFrags.fastq -out ' + path.replace('overlapped/', 'trimmed/') + 
           sample + ' -log ' + path.replace('overlapped/','') + 'logs/primer_trimming/' + 
           sample + '.primertrim.log -nomatch 3 -tag5 ACTCCTACGGGAGGCAGCAG -mm5 2' + 
           ' -tag3 ATTAGAWACCCBDGTAGTCC -mm3 2')
    run(cmd)

# Check for empty files
def check_files(path, f, tag):
    n = os.path.getsize(path + f + tag)
    return n

# Quality filter with usearch8 
def qc(path, f):
    cmd = ('../tools/usearch8 -fastq_filter ' + path + f + '.fastq ' + 
           '-fastaout ' + path.replace('/trimmed/', '/quality_filtered/') + f +
           '.fasta -fastq_maxee 6 -relabel ' + f + ' -eeout -quiet')
    run(cmd)

# Generate a make.contigs() mothur command with all samples
def gen_makecontigs(path):
    os.chdir(path)
    files = glob.glob('*.fasta')
    mothurcmd = 'make.groups(fasta='
    for f in files:
        mothurcmd = mothurcmd + f + '-'
    mothurcmd = mothurcmd.rstrip('-') + ', groups='
    for f in files:
        mothurcmd = mothurcmd + f.replace('.fasta', '') + '-'
    mothurcmd = mothurcmd.rstrip('-') + ')'
    return mothurcmd

if __name__ == '__main__':
  main()
