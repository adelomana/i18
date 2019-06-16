#!/bin/bash
freebayes --fasta-reference /Volumes/omics4tb2/alomana/software/cellRanger/data/refdata-cellranger-mm10-3.0.0/fasta/genome.fa /Volumes/omics4tb2/alomana/projects/i18/results/deconvolution/mouse/d4_pstat_r2_mouse/outs/possorted_genome_bam.bam --min-alternate-count 10 --min-alternate-fraction 1.0 --min-coverage 10 -g 200 --min-mapping-quality 30 --min-base-quality 20 > freebayes.d4_pstat_r2_mouse.messages.vcf
