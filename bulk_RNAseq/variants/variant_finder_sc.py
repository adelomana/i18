import os,datetime,sys,pickle
import multiprocessing,multiprocessing.pool

def cell_analyzer(argument):

    cellID=argument[0]
    cellIDs=argument[1]

    printt('working with cell {}/{}'.format(cellIDs.index(cellID),len(cellIDs)))

    # f.1. determine the bam file to inspect given the cellID (which well of the 10x)
    if cellID[-1] == '1':
        tag='Time0_1-WT_1-MAVS_1-NLRP3'
    elif cellID[-1] == '2':
        tag='Time4_1-WT_1-MAVS_1-NLRP3'
    elif cellID[-1] == '3':
        tag='Time24_1-WT_1-MAVS_1-NLRP3'
    elif cellID[-1] == '4':
        tag='Time24b_1-WT_1-MAVS_6-NLRP3'
    else:
        raise ValueError('well ID not found')
    
    bam_file=bam_dir+tag+'/outs/possorted_genome_bam.bam'

    # f.2. define the output dir
    main_dir=variants_dir+tag+'/'
    if os.path.exists(main_dir) == False:
        os.mkdir(main_dir)

    working_dir=main_dir+'{}/'.format(cellID)
    if os.path.exists(working_dir) == False:
        os.mkdir(working_dir)
    
    sc_bam_file=working_dir+cellID+'.bam'

    # f.3. defined barcode file
    one_cell_barcode_file=working_dir+'one_cell_barcode.txt'
    with open(one_cell_barcode_file,'w') as f:
        f.write('{}\n'.format(cellID))
    f.close()

    # f.3. build the command    
    flag1='subset-bam'
    flag2='--bam {}'.format(bam_file)
    flag3='--cell-barcodes {}'.format(one_cell_barcode_file)
    flag4='--out-bam {}'.format(sc_bam_file)
    flag5='--cores 8'

    flags=[flag1,flag2,flag3,flag4]
    command=' '.join(flags)

    # f.4. run the command    
    printt(command)
    os.system(command)

    # f.5. check number of mapped reads
    mapped_reads_file=working_dir+'mapped_reads.txt'
    command='samtools view -c -F 260 {} > {}'.format(sc_bam_file,mapped_reads_file)
    printt(command)
    os.system(command)

    # f.6. call variants
    target_vcf_file=working_dir+'/one_cell_variants.vcf'
    command='freebayes --fasta-reference {} {} --min-alternate-count 2 --min-alternate-fraction 0.66 --skip-coverage 200 --min-mapping-quality 30 --min-base-quality 20 > {}'.format(fasta_reference_file,sc_bam_file,target_vcf_file)
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
selected_cellIDs_file='/proj/omics4tb2/alomana/projects/i18/results/10x/selectedIDs/selected_ids.pickle'
bam_dir='/proj/omics4tb2/alomana/projects/i18/results/10x/'
fasta_reference_file='/proj/omics4tb2/alomana/software/cellRanger/data/refdata-cellranger-mm10-3.0.0/fasta/genome.fa'
variants_dir='/proj/omics4tb2/alomana/projects/i18/results/10x/variants/'
threads=64

# 1. recover selected cellIDs
printt('recovering cellIDs')
f=open(selected_cellIDs_file,'rb')
cellIDs=pickle.load(f)
f.close()
printt('recovered {} selected cell IDs'.format(len(cellIDs)))

# 2. process single cells
arguments=[[cellID,cellIDs] for cellID in cellIDs]
hydra=multiprocessing.pool.Pool(threads)
empty=hydra.map(cell_analyzer,arguments)
hydra.close()

"""
for cellID in cellIDs[:8]:
    argument=[cellID,cellIDs]
    cell_analyzer(argument)
""" 
