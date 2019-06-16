###
### This script will generate cell-specific bam files and will store them in separate directories sorted by number of mapped number of reads
###

import os,sys
import multiprocessing,multiprocessing.pool

def barcode_processor(barcode):

    print(barcode)

    # create a directory with a barcode
    target_dir=output_dir+barcode+'/'
    target_cell_barcode_file=target_dir+barcode+'.tsv'
    target_bam_file=target_dir+barcode+'.bam'
    
    if os.path.exists(target_dir) == False:
        os.mkdir(target_dir)
        with open(target_cell_barcode_file,'w') as f:
            f.write(barcode)

    # build the subset-bam command
    command='subset-bam --bam {} --cell-barcodes {} --out-bam {}'.format(chromium_bam_location,target_cell_barcode_file,target_bam_file)
    ###print('\t',command)
    os.system(command)

    # define the number of mapped reads
    mapped_reads_file=target_dir+'mapped_reads.txt'
    command='samtools view -c -F 260 {} > {}'.format(target_bam_file,mapped_reads_file)
    #print('\t',command)
    os.system(command)

    return None

###
### MAIN
###

# 0. user defined variables
barcodes_file='/Users/alomana/scratch/d4_pstat_r2_mouse/outs/filtered_feature_bc_matrix/barcodes.tsv'
chromium_bam_location='/Users/alomana/scratch/d4_pstat_r2_mouse/outs/possorted_genome_bam.bam'
output_dir='/Users/alomana/scratch/d4_pstat_r2_mouse/outs/cell_bam_files/'
threads=8

# 1. generate cell-specific BAM files

# 1.1. read barcodes
barcodes=[]
with open(barcodes_file,'r') as f:
    for line in f:
        v=line.split()[0]
        barcodes.append(v)

# 1.2. generate bam files
#barcodes=barcodes[:5]       
#for barcode in barcodes:
#    barcode_processor(barcode)

hydra=multiprocessing.pool.Pool(threads)
labels=hydra.map(barcode_processor,barcodes)
hydra.close()

# 2. identify which cells have more aligned reads
