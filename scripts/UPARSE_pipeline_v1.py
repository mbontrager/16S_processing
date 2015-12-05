#!/usr/bin/python

import os, sys, getopt, glob, subprocess, shutil

############################################################
# Run the UPARSE pipeline on demultiplexed, overlapped, 
# fasta files. 
#
# Full path to "quality_filtered" directory from demux_samples.py
# Requires UPARSE gold database, greengenes fasta, greengenes taxonomy in 
# "../database" (relative to the "quality_filtered" directory). 
# These are relatively large databases (GG especially).
# email mbontrager@gmail.com if you need the databases.
#
# Usage: Run from script directory:
# UPARSE_pipeline.py -p path/to/fasta/dir
#
# Author: Martin Bontrager
############################################################

def main():
    path = ''
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],'hp:',['path='])
    except getopt.GetoptError:
        print('UPARSE_pipeline.py -p path/to/fasta/dir')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('UPARSE_pipeline.py -p path/to/fasta/dir')
            sys.exit()
        elif opt in ('-p', '--path'):
            path = arg
    
    os.chdir(path)
    sample_list = glob.glob('*.fasta')
    add_name_to_header(sample_list)
    run('mkdir fixed_headers')
    run('mv *.fa fixed_headers')

    # Pool samples into one large file
    os.chdir('fixed_headers')
    outfilename = 'concatenated.fa'
    with open(outfilename, 'wb') as outfile:
        for filename in glob.glob('*.fa'):
            if filename == outfilename:
                continue
            with open(filename, 'rb') as readfile:
                shutil.copyfileobj(readfile, outfile)
    run('mv ' + outfilename + ' ..')
    os.chdir('..')
    run('rm -r fixed_headers/')
    
    # Screen for short/long sequences:
    cmd = ('mothur "#screen.seqs(fasta=concatenated.fa,' + 
           ' minlength=390, maxlength=450, processors=5)"')
    run(cmd)
    run('rm concatenated.bad.accnos concatenated.fa')
    # Dereplicate sequences and sort by binned size, discard singletons(
    # See UPARSE documenation)
    run('usearch -derep_fulllength concatenated.good.fa -fastaout' + 
        ' uniques.fasta -sizeout -threads 2 -relabel BIN')
    run('usearch -sortbysize uniques.fasta -fastaout seqs_sorted.fasta' +
        ' -minsize 2')
    # Cluster into OTUs @ 97% ID
    run('usearch -cluster_otus seqs_sorted.fasta -otus otus.fa' + 
        ' -uparseout results.txt -relabel OTU_ -sizeout')
    # Filter out remaining chimeric reads and generate output table
    run('usearch -uchime_ref otus.fa -db ../database/uchime_gold.fa ' + 
        '-nonchimeras otu.uchime.fa -strand plus -threads 2')
    # Perform a series of steps to classify OTUs by taxonomy with UTAX
    cmd = ('usearch -utax otu.uchime.fa -db ' +
           '../database/tax.udb -utaxout reads.utax ' +
           '-strand both -alnout aln.txt')
    run(cmd)
    run('usearch -usearch_global concatenated.good.fa -db otu.uchime.fa ' + 
        '-uc readmap.uc -strand both -id 0.97')
    run('python /home/lee/bioinformatics/drive5/uc2otutab.py readmap.uc' + 
        ' > otutable.txt')
    run('rm uniques.fasta seqs_sorted.fasta otus.fa concatenated.good.fa')
    run('mothur "#align.seqs(fasta=otu.uchime.fa, ' + 
        'reference=../database/silva.v3v4.fasta)"')
    run('mothur "#screen.seqs(fasta=otu.uchime.align, criteria=90,' + 
        ' optimize=start-end)"')
    run('mothur "#filter.seqs(fasta=otu.uchime.good.align,' + 
        ' vertical=T, trump=.)"')
    run('mv otu.uchime.good.filter.fasta alignment.fasta')
    run('/media/DATA/programs/FastTreeMP -gtr -nt alignment.fasta' + 
        ' > FastTree2.tre')
    run('sed "s/\(OTU_[0-9]*\);size=[0-9]*;/\\1/g" <FastTree2.tre' + 
        ' >ParsedFastTree2.tre')
    # Clean and organize
    run('mkdir ../2015-11-17_UPARSE')
    run('mv otutable.txt alignment.fasta ParsedFastTree2.tre otus_tax.fa ' +
        'aln.txt FastTree2.tre ../2015-11-17_UPARSE')
    run('mv results.txt ../2015-11-17_UPARSE')
    run('rm mothur* otu* readmap.uc')
 
def add_name_to_header(files):
    for i in files:
        name = i.replace('.fasta', '')
        cmd = ('sed \"-es/^>\(.*\)/>\\1;barcodelabel=' + name + ';/\" < ' + 
        i + ' > ' + name + '.fa')
 
        run(cmd)
        
# Simplify running bash commands
def run(cmd):
    print('Running command: \n' + cmd + '\n\n') 
    p = subprocess.Popen(cmd, shell=True)
    p.communicate()
      
if __name__ == "__main__":
    main()
