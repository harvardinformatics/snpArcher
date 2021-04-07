#!/bin/bash
#SBATCH -J sm
#SBATCH -o out
#SBATCH -e err
#SBATCH -p holy-info
#SBATCH -n 1
#SBATCH -t 9000
#SBATCH --mem=10000


source activate snakemake
snakemake --snakefile Snakefile_fastq2bam --profile ./profiles/slurm
