#!/bin/bash
#SBATCH -J sm
#SBATCH -o out
#SBATCH -e err
#SBATCH -p 256x20
#SBATCH -n 1
#SBATCH -t 9000
#SBATCH --mem=10000

snakemake --snakefile workflow/Snakefile --profile profiles/slurm
