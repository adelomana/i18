#!/bin/bash

#$ -N sc_finder
#$ -o /proj/omics4tb2/alomana/scratch/messages.sc.calling.o.txt
#$ -e /proj/omics4tb2/alomana/scratch/messages.sc.calling.e.txt
#$ -pe smp 64
#$ -S /bin/bash

cd /users/alomana
source .bash_profile 

echo "starting..."
cd /proj/omics4tb2/alomana/projects/i18/deployment
time python variant_finder_sc.py > messages.txt
echo "... run completed."
