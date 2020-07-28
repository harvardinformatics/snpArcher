#!/bin/bash
#SBATCH -J sm
#SBATCH -o out
#SBATCH -e err
#SBATCH -p holy-info
#SBATCH -n 1
#SBATCH -t 4000
#SBATCH --mem=4000


#snakemake --snakefile Snakefile_bam2vcf --profile ./profiles/slurm
snakemake --snakefile Snakefile_bam2vcf_fb --profile ./profiles/slurm

