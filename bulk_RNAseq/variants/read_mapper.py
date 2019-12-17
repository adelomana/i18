###
### this script runs the mapping of reads
###

import os,sys

def genome_index_maker():

    flag1=' --runMode genomeGenerate'
    flag2=' --runThreadN %s'%threads
    flag3=' --genomeDir %s'%genome_index_dir
    flag4=' --genomeFastaFiles %s'%genome_fasta
    flag5=' --sjdbGTFfile %s'%genome_annotation
    flag6=' --sjdbOverhang 74'

    cmd=STAR_executable+flag1+flag2+flag3+flag4+flag5+flag6

    print()
    print(cmd)
    print()
    os.system(cmd)

    return None

def STAR_calling(genotype):

    read_files=[]
    for label in genotypes[genotype]:
        for i in range(4):
            read_files.append(reads_dir+label+'_L00'+str(i+1)+'.clean.fastq')
    read_paths=','.join(read_files)

    out_dir=bam_dir+genotype+'/'
    if os.path.exists(out_dir) == False:
        os.mkdir(out_dir)

    flag1=' --genomeDir %s'%genome_index_dir
    flag2=' --runThreadN %s'%threads
    flag3=' --readFilesIn %s'%read_paths 
    flag4=' --outFileNamePrefix {}'.format(out_dir)
    flag5=' --outSAMtype BAM SortedByCoordinate'
    flag6=' --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000' # ENCODE options

    cmd='time '+STAR_executable+flag1+flag2+flag3+flag4+flag5+flag6
    
    print()
    print(cmd)
    print()
    os.system(cmd)
    
    return None

###
### MAIN
###

# 0. user defined variables
reads_dir='/proj/omics4tb2/alomana/projects/i18/results/bulk_rnaseq/clean_fastq/'
bam_dir='/proj/omics4tb2/alomana/projects/i18/results/bulk_rnaseq/bam/'
STAR_executable='/proj/omics4tb2/alomana/software/STAR-2.7.3a/bin/Linux_x86_64/STAR'
threads=16
genome_index_dir='/proj/omics4tb2/alomana/projects/i18/results/bulk_rnaseq/star_index'
genome_fasta='/proj/omics4tb2/alomana/software/cellRanger/data/refdata-cellranger-mm10-3.0.0/fasta/genome.fa'             
genome_annotation='/proj/omics4tb2/alomana/software/cellRanger/data/refdata-cellranger-mm10-3.0.0/genes/genes.gtf'

# 1. generate index
#print('generate index...')
#genome_index_maker()

# 2. retreieve samples
files=os.listdir(reads_dir)
labels=[element.split('_L')[0] for element in files]
labels=list(set(labels))
labels.sort()

genotypes={}
genotypes['wt']=[]; genotypes['nlrp3']=[]; genotypes['mavs']=[]
for label in labels:
    if 'wt' in label:
        genotypes['wt'].append(label)
    if 'NLRP3' in label:
        genotypes['nlrp3'].append(label)
    if 'MAVS' in label:
        genotypes['mavs'].append(label)

# 3. map reads
for genotype in genotypes:
    STAR_calling(genotype)
