import os,datetime,sys
import multiprocessing,multiprocessing.pool

def printt(message):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))

    return None

###
### MAIN
###

# 0. user defined variables
selected_cellIDs_file='xx'
bam_dir='/Volumes/omics4tb2/alomana/projects/i18/results/bulk_rnaseq/bam/'
fasta_reference_file='/Volumes/omics4tb2/alomana/software/cellRanger/data/refdata-cellranger-mm10-3.0.0/fasta/genome.fa'
threads=3

