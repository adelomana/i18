#!/bin/bash

#$ -N mapping
#$ -o /proj/omics4tb2/alomana/scratch/messages.mapping.o.txt
#$ -e /proj/omics4tb2/alomana/scratch/messages.mapping.e.txt
#$ -pe smp 16
#$ -l hostname=baliga1
#$ -S /bin/bash

cd /users/alomana
source .bash_profile 

echo "starting..."
cd /proj/omics4tb2/alomana/projects/i18/src/bulk
time python read_mapper.py
echo "... run completed."
