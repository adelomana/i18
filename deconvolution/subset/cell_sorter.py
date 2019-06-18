###
### this script sorts cells based on snps
### provides output as genotype, cell id, number of mapped reads, total variants found, number of reference snps, number of genotype specific snps (3), proportions of genotype specific snps (3), prediction 
###

import sys,datetime
import multiprocessing,multiprocessing.pool

def cell_classifier(cellID):

    # read the number of mapped reads
    info_file=chromium_folder+'outs/cell_bam_files/'+cellID+'/mapped_reads.txt'
    with open(info_file,'r') as f:
        mapped_reads=int(f.readline())

    # read the variants
    detected_variants=0
    detected_snps=0
    condition_snp_counts=[0 for condition in conditions] # alphabetically sorted  

    variants_file=chromium_folder+'outs/cell_bam_files/'+cellID+'/'+cellID+'.vcf'
    try:
        with open(variants_file,'r') as f:
            for line in f:
                if line[0] != '#':
                    detected_variants=detected_variants+1
                    v=line.split('\t')
                    if len(v[4]) == 1:
                        detected_snps=detected_snps+1
                        label='{}.{}.{}'.format(v[0],v[1],v[4])

                        for i in range(len(conditions)):
                            if label in sorted_snps[conditions[i]]:
                                condition_snp_counts[i]=condition_snp_counts[i]+1

                                message='{} {} {} {} {}'.format(cellID,mapped_reads,condition_snp_counts,conditions[i],label)
                                printt(message)
    except:
        message='{} {} vcf file not found'.format(cellID,mapped_reads)
        printt(message)


    print(cellID,mapped_reads,detected_snps)

                
                

    return None

def printt(message):
    
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))

    return None

###
### MAIN
###

# 0. user defined variables
threads=8
genotypes=['wt','mavs','nlrp3']
snps_dir='/Volumes/omics4tb2/alomana/projects/i18/results/deconvolution/freebayes/'

working_genotype=genotypes[0]
chromium_folder='/Volumes/omics4tb2/alomana/projects/i18/results/deconvolution/mouse/d2_dmso_mouse/'

printt('working with {} genotype'.format(working_genotype))

# 1. read reference and genotype
printt('read pre-defined snps')
sorted_snps={}
conditions=[element for types in [genotypes,['reference']] for element in types]
conditions.sort()

for condition in conditions:
    sorted_snps[condition]=[]
    vcf_file=snps_dir+condition+'.vcf'
    with open(vcf_file,'r') as f:
        for line in f:
            v=line.split()[0]
            sorted_snps[condition].append(v)

# 2. sort cells

# 2.1. define cell barcodes
barcode_file=chromium_folder+'outs/filtered_feature_bc_matrix/barcodes.tsv'
cellIDs=[]
with open(barcode_file,'r') as f:
    for line in f:
        v=line.split()[0]
        cellIDs.append(v)
printt('defined {} cell barcodes'.format(len(cellIDs)))

# 2.2. retrieve cell-specific variants
printt('sort cells')

#cellIDs=cellIDs[:8]

hydra=multiprocessing.pool.Pool(threads)
labels=hydra.map(cell_classifier,cellIDs)
hydra.close()

#for cellID in cellIDs:
#    try:
#        cell_classifier(cellID)
#    except:
#        print('{} \t did not work'.format(cellID))
