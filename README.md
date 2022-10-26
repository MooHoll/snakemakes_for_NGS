# Snakemake workflows on ALICE2: Univeristy of Leicester HPC

Use a conda environment to call snakemake, the module doesn't seem to be set up properly. Create a config file, Snakefile and .PBS script.


conda create --name snakemake_env

conda activate snakemake_env

conda install -c bioconda snakemake

NOTE: If running tests need to run: snakemake --unlock when done, before submitting the job to the cluster,
otherwise it will kill out.