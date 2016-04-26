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

# Change to True if input files are larger than 4 GB (limit for usearch)
subsample = False

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

    # Sub-sample fasta files since the latest batch is too big to process
    if subsample:
        for f in sample_list:
            cmd = ('mothur "#sub.sample(fasta=' + f + ')"')
            run(cmd)
        sample_list = glob.glob('*.fasta')
    add_name_to_header(sample_list)
    run('mkdir fixed_headers')
    run('mv *.fa fixed_headers')

    pool()
    screen()
    usearch()
    mothur()

def add_name_to_header(files):
    for i in files:
        name = i.replace('.subsample.fasta', '')
        name = name.replace('.fasta', '')
        cmd = ('sed \"-es/^>\(.*\)/>\\1;barcodelabel=' + name + ';/\" < ' + 
        i + ' > ' + name + '.fa')
 
        run(cmd)    

# Pool sequences and process
def pool():
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
def screen():
    cmd = ('mothur "#screen.seqs(fasta=concatenated.fa,' + 
           ' minlength=390, maxlength=450, processors=5)"')
    run(cmd)
    run('rm concatenated.bad.accnos concatenated.fa')
        
# Mothur taxonomic classification and tree building
def mothur():
    ''' Perform a series of steps to classify OTUs by taxonomy with Greengenes
    and align seqs in multiple alignment to create FastTree tree.'''
    
    cmd = ('mothur "#classify.seqs(fasta=otu.uchime.fa, ' + 
          'template=../database/gg_13_8_FWSET.fasta, processors=5, ' + 
          'taxonomy=../database/gg_13_8_FWSET.tax, cutoff=70)"')
    run(cmd)
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
    run('mkdir ../UPARSE')
    run('mv otutable.txt alignment.fasta ParsedFastTree2.tre ../UPARSE')
    run('mv FastTree2.tre ../UPARSE')
    run('mv otu.uchime.gg* results.txt ../UPARSE')
    run('rm mothur* otu* readmap.uc')

# Run USEARCH/UCLUST
def usearch():
    # Dereplicate sequences and sort by binned size
    run('usearch -derep_fulllength concatenated.good.fa -fastaout' + 
        ' uniques.fasta -sizeout -threads 2 -relabel BIN')
    run('usearch -sortbysize uniques.fasta -fastaout seqs_sorted.fasta')
    # Cluster into OTUs @ 97% ID
    run('usearch -cluster_otus seqs_sorted.fasta -otus otus.fa' + 
        ' -uparseout results.txt -relabel OTU_ -sizeout')
    # Filter out remaining chimeric reads and generate output table
    run('usearch -uchime_ref otus.fa -db ../database/uchime_gold.fa ' + 
        '-nonchimeras otu.uchime.fa -strand plus -threads 2')
    run('usearch -usearch_global concatenated.good.fa -db otu.uchime.fa ' + 
        '-strand plus -id 0.97 -uc readmap.uc -maxaccepts 8 -maxrejects 64 ' + 
        '-top_hit_only')
    run('rm uniques.fasta seqs_sorted.fasta otus.fa concatenated.good.fa')
    run('python /home/lee/bioinformatics/drive5/uc2otutab.py readmap.uc' + 
        ' > otutable.txt')
        
# Simplify running bash commands
def run(cmd):
    print('Running command: \n' + cmd + '\n\n') 
    p = subprocess.Popen(cmd, shell=True)
    p.communicate()
    
if __name__ == "__main__":
    main()
