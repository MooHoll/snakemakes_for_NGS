#!/bin/bash

#PBS -N snakemake
#PBS -l walltime=24:00:00
#PBS -l vmem=50gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=20

# Run script in the working directory it was submitted in  (8 nodes and 6hrs)
cd $PBS_O_WORKDIR 

source ~/miniconda3/bin/activate snakemake_env

snakemake