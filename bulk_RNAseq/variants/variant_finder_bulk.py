import os,datetime,sys
import multiprocessing,multiprocessing.pool

def freebayes_caller(label):

    bam_file='{}{}/Aligned.sortedByCoord.out.bam'.format(bam_dir,label)
    target_vcf_file='{}.vcf'.format(label)

    command='freebayes --fasta-reference {} {} --min-alternate-count 2 --min-mapping-quality 30 --min-base-quality 20 > {}'.format(fasta_reference_file,bam_file,target_vcf_file)

    printt(command)
    os.system(command)

    return None

def printt(message):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))

    return None

###
### MAIN
###

# 0. user defined variables
bam_dir='/Volumes/omics4tb2/alomana/projects/i18/results/bulk_rnaseq/bam/'
fasta_reference_file='/Volumes/omics4tb2/alomana/software/cellRanger/data/refdata-cellranger-mm10-3.0.0/fasta/genome.fa'
threads=3

# 1. define the sets of BAM files per genotype
folders=os.listdir(bam_dir)

# 2. call freebayes

#for folder in folders:
#    freebayes_caller(folder)

hydra=multiprocessing.pool.Pool(threads)
empty=hydra.map(freebayes_caller,folders)
hydra.close()
