#!/bin/bash
#SBATCH -J sm
#SBATCH -o out
#SBATCH -e err
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 9000
#SBATCH --mem=10000


source activate snakemake
snakemake --snakefile Snakefile_intervals --profile ./profiles/slurm

snakemake --snakefile Snakefile_bam2vcf_gatk --profile ./profiles/slurm --dryrun > bam2vcf_gatk_dryrun.txt
snakemake --snakefile Snakefile_bam2vcf_fb --profile ./profiles/slurm --dryrun > bam2vcf_fb_dryrun.txt
