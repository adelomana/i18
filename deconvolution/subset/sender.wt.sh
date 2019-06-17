#!/bin/bash

#$ -N wt
#$ -o /proj/omics4tb2/alomana/scratch/messages.wt.o.txt
#$ -e /proj/omics4tb2/alomana/scratch/messages.wt.e.txt
#$ -pe smp 64
#$ -S /bin/bash

cd /users/alomana
source .bash_profile 

echo "Original paths given my bash_profile..."
echo $PATH
echo ""

echo "running pipeline..." 
cd /proj/omics4tb2/alomana/projects/i18/src
python pipeline.v2.wt.py
echo "... run completed"
