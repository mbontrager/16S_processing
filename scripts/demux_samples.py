#!/usr/bin/python

import os, sys, getopt, subprocess, glob

############################################################
# Demultiplex MiSeq Reads with internal barcode
#
# Pipeline: a QIIME combined seqs file by sample. Requires
# local QIIME installation.
#
# Usage: demux_samples.py -1 read1.fastq -2 read2.fastq
#
# Author: Martin Bontrager
############################################################

def main():
    read1 = ''
    read2 = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"h1:2:",["r1=","r2="])
    except getopt.GetoptError:
        print('demux_samples.py -1 <read1.fastq> -2 <read2.fastq>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('demux_samples.py -1 <read1.fastq> -2 <read2.fastq>')
            sys.exit()
        elif opt in ("-1", "--r1"):
            read1 = arg
        elif opt in ("-2", "--r2"):
            read2 = arg
    path1 = os.path.dirname(read1) + '/'
    path2 = os.path.dirname(read2) + '/'
    get_barcodes(read1, read2, path1, path2)
    trim_barcodes(read1, read2, path1, path2)
    
# Simplify running bash commands
def run(cmd):
    p = subprocess.Popen(cmd, shell=True)
    os.waitpid(p.pid, 0)

# Get barcodes from both reads
def get_barcodes(read1, read2, path1, path2):
    cmd = ('../tools/fastx_trimmer -i ' + 
            read1 + ' -f 1 - l 12 -Q 33 -o ' + path1 + 'R1_barcode.fq')
    run(cmd)
    cmd = ('../tools/fastx_trimmer -i ' + 
            read2 + ' -f 1 - l 12 -Q 33 -o ' + path2 + 'R2_barcode.fq')
    run(cmd)
    cmd = ('cat ' + path1 + 'R1_barcode.fq | '
            '../tools/MiSeq16S-master/fq_mergelines.pl >'
             + path1 + 'R1_barcode_temp')
    run(cmd)
    cmd = ('cat ' + path2 + 'R2_barcode.fq | '
            '../tools/MiSeq16S-master/fq_mergelines.pl >'
            + path2 + 'R2_barcode_temp')
    run(cmd)
    cmd = ('paste ' + path1 + 'R1_barcode_temp ' + path2 + 'R2_barcode_temp |'
           ' awk -F"\\t" '
          """'{print $5"\\t"$2$6"\\t"$3"\\t"$4$8}' """
          "| ../tools/MiSeq16S-master/fq_splitlines.pl > " +
          path1 + "R1R2_barcode.fastq")
    run(cmd)

# Trim barcodes from original sequence files
def trim_barcodes(read1, read2, path1, path2):
    cmd = ('../tools/seqtk/seqtk trimfq -b 12 ' + read1 + ' > ' + 
            path1 + 'R1_trimmed_seq.fastq')
    run(cmd)
    cmd = ('../tools/seqtk/seqtk trimfq -b 12 ' + read2 + ' > ' + 
            path2 + 'R2_trimmed_seq.fastq')
    run(cmd)


if __name__ == "__main__":
    main()