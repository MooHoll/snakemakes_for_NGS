# snakemakes_for_NGS

On ALICE2 HPC in Leicester.

conda create --name snakemake_env
conda activate snakemake_env
conda install -c bioconda snakemake

snakemake --cores 2 --use_conda 

