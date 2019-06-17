###
### this script finds the variants specific of each genotype
###

import os,sys,datetime

def printt(message):
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))
    return None

# 0. user defined variables
genotypes_dir='/Volumes/omics4tb2/alomana/projects/i18/results/deconvolution/freebayes/'
genotypes_files={}
genotypes_files['wt']='freebayes.d2_dmso_mouse.messages.vcf'
genotypes_files['nlrp3']='freebayes.d4_pstat_r2_mouse.messages.vcf'
genotypes_files['mavs']='freebayes.d4_dmso_mouse.messages.vcf'
genotypes=list(genotypes_files.keys())
genotypes.sort()

# 1. read variants per genotype
variants={} # variants[genotype]=['chr.pos.var',original line]
all_snps_lists=[]

for genotype in genotypes:

    printt('working on genotype {}'.format(genotype))
    
    detected_variants=0
    detected_snps=0
    all_snps_lists.append([])
    
    variants[genotype]=[]

    vcf_file=genotypes_dir+genotypes_files[genotype]
    with open(vcf_file,'r') as f:
        for line in f:
            if line[0] != '#':
                detected_variants=detected_variants+1
                v=line.split('\t')
                if len(v[4]) == 1:
                    detected_snps=detected_snps+1
                    label='{}.{}.{}'.format(v[0],v[1],v[4])
                    if label not in all_snps_lists[-1]:
                        all_snps_lists[-1].append(label)
                    variants[genotype]=[label,line]
    printt('detected {} variants'.format(detected_variants))
    printt('detected {} snps'.format(detected_snps))
    
# 2. define specific and reference snps
printt('join genotype snps')
all_snps=[]
for elements in all_snps_lists:
    for snp in elements:
        if snp not in all_snps:
            all_snps.append(snp)
printt('working with {} total snps'.format(len(all_snps)))

sorted_snps={}
for genotype in genotypes:
    sorted_snps[genotype]=[]
number_reference_snps=0
    
for snp in all_snps:
    found=[]
    for i in range(len(all_snps_lists)):
        if snp in all_snps_lists[i]:
            found.append(i)
    if len(found) == 1:
        sorted_snps[genotypes[found[0]]].append(snp)
    if len(found) == 3:
        number_reference_snps=number_reference_snps+1

for genotype in genotypes:
    printt('found {} especific snps for {}'.format(len(sorted_snps[genotype]),genotype))
printt('found {} reference snps'.format(number_reference_snps))

# 3. store genotype specific variants
printt('write genotype especific snps')
for genotype in sorted_snps:
    storage_file=genotype+'.vcf'
    with open(storage_file,'w') as f:
        for snp in sorted_snps[genotype]:
            f.write('{}\n'.format(snp))
