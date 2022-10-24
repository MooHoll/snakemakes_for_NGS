#!/bin/bash

#PBS -N getting_data
#PBS -l walltime=00:05:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8
#PBS -q devel

# Change directory to the one the job was submitted in
cd $PBS_O_WORKDIR 

# Load required modules
module load sratoolkit/2.11.1

#prefetch --option-file accessions.txt

#move all files out of the directories

cat accessions.txt | while read line
do
fasterq-dump --split-files ${line}.sra
done