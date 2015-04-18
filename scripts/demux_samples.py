#!/usr/bin/python

import os, sys, getopt, subprocess, glob, re
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
    proc = '1'
    scriptdir = os.getcwd()

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hkf:r:c:p:",["forward=","reverse=", "csvfile=", "processors="])
    except getopt.GetoptError:
        print('demux_samples.py -f <forward_read.fastq> -r <reverse_read.fastq> -c <map.csv>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('demux_samples.py -f <forward_read.fastq> -r <reverse_read.fastq> -c <map.csv>')
            sys.exit()
        if opt == '-k':
            keep = True
        elif opt in ("-p", "--processors"):
            proc = str(arg)
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
    subprocess.call(['mkdir', (path + 'mothur')])

    change_names(csvfile)
    csvfile = csvfile.replace('.csv', '_fixed.csv')
    csv_parse(csvfile, (path + 'barcodes.fil'), (path + 'samples.txt'))
    demux(forw, rev, path)
    
    for i in get_samples(path + 'samples.txt'):
        if check_files(dpath, i, '_338F.fastq') == 0:
            continue        
        trim_barcodes((dpath + i + '_338F.fastq'), (dpath + i + '_806R.fastq'))
        overlap_usearch(i, (dpath + i + '_338F_bctrimmed.fastq'), (dpath + i + '_806R_bctrimmed.fastq'))
    print('Overlap complete')
    clean(keep, dpath)
    
    for i in get_samples(path + 'samples.txt'):
        if not os.path.isfile(path + 'overlapped/' + i + '.extendedFrags.fastq'):
            continue     
        elif check_files((path + 'overlapped/'), i, '.extendedFrags.fastq') == 0:
            continue
        trim_primers((path + 'overlapped/'), i)
    print('Primers trimmed')
    clean(keep, (path + 'overlapped/'))
    
    for i in glob.glob(path + 'trimmed/*.fastq'):
        qc(i)
    print('QC complete. WARNING errors generally indicate one or more empty samples')
    clean(keep, (path + 'trimmed/'))
    mothur(path, proc, scriptdir)

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

# Overlap reads using flash (NOT USED)
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
           reverse + ' -fastqout '+ d + sample + '.extendedFrags.fastq -log ' + 
           d.replace('/overlapped/', '/logs/overlap/') + sample + '.overlap.log -quiet')
    run(cmd)

# Trim 338F/806R primers with tagcleaner
def trim_primers(path, sample):
    cmd = ('../tools/tagcleaner-standalone-0.16/tagcleaner.pl -fastq ' + path + 
           sample + '.extendedFrags.fastq -out ' + path.replace('overlapped/', 'trimmed/') + 
           sample + ' -nomatch 3 -tag5 ACTCCTACGGGAGGCAGCAG -mm5 2' + 
           ' -tag3 ATTAGAWACCCBDGTAGTCC -mm3 2')
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    p.communicate()


# Check for empty files
def check_files(path, f, tag):
    n = os.path.getsize(path + f + tag)
    return n

# Quality filter with usearch8 
def qc(f):
    cmd = ('../tools/usearch8 -fastq_filter ' + f + 
           ' -fastaout ' + f.replace('/trimmed/', '/quality_filtered/').replace('.fastq', '.fasta') + 
           ' -fastq_maxee 6 -quiet')
    run(cmd)

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
def mothur(path, p, s):
    f = open('batch.mothur', 'r')
    tmp = f.readlines()
    f.close()
    f = open((path + 'batch.mothur'), 'w')
    for i, line in enumerate(tmp):
        f.write(line)
        if i == 0:
            f.write('set.dir(input=' + path + 'quality_filtered)\n')
            f.write(gen_makegroups(path + 'quality_filtered'))
            f.write('set.dir(input=' + path + 'mothur, output='+ path + 'mothur)\n')
            f.write('system(mv ' + path + 'quality_filtered/mergegroups ' + path + 'mothur)\n')
            f.write('summary.seqs(fasta=allsamples.fasta, processors=' + p + ')\n')
    f.close()

    os.chdir(path + 'mothur/')
    run('mothur ../batch.mothur')

# Clean up temp files from demultiplexing and trimming
def clean(k, f):
    if k is False:
        subprocess.call(['rm', '-r', f])
    else:
        pass

if __name__ == '__main__':
  main()
