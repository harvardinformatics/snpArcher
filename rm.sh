#!/bin/bash
#SBATCH -J rm
#SBATCH -o out
#SBATCH -e err
#SBATCH -p holy-info
#SBATCH -n 1
#SBATCH -t 40
#SBATCH --mem=1000

rm -r 00_fastqFiltered
rm -r 01_mappedReads
rm -r 02_bamSumstats
rm -r 03_gvcfs
rm -r 04_genomicsDB
rm -r 05_vcfs
rm Combined*.vcf
rm Combined*.idx
rm slurm-*.out
rm out
rm err
rm bam_sumstats.txt
