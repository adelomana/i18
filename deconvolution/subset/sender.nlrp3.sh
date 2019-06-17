#!/bin/bash

#$ -N nlrp3
#$ -o /proj/omics4tb2/alomana/scratch/messages.nlrp3.o.txt
#$ -e /proj/omics4tb2/alomana/scratch/messages.nlrp3.e.txt
#$ -pe smp 40
#$ -S /bin/bash

cd /users/alomana
source .bash_profile 

echo "Original paths given my bash_profile..."
echo $PATH
echo ""

echo "running pipeline..." 
cd /proj/omics4tb2/alomana/projects/i18/src
python pipeline.v2.nlrp3.py
echo "... run completed"
