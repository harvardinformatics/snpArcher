#!/bin/bash
#SBATCH -J rm
#SBATCH -o out
#SBATCH -e err
#SBATCH -p holy-info
#SBATCH -n 1
#SBATCH -t 40
#SBATCH --mem=1000

rm -r 01_mappedReads
rm -r 02_sumstats
rm slurm-*.out
rm out
rm err
rm sumstats.txt
