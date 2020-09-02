#!/bin/bash
#SBATCH -J rm
#SBATCH -o out
#SBATCH -e err
#SBATCH -p holy-info
#SBATCH -n 1
#SBATCH -t 40
#SBATCH --mem=1000

rm -r fastq2bam
rm -r intervalFiles
rm -r gatk
rm -r freebayes

rm -r logs

rm data/zebraFinch/genome/Tgut_subseg_renamed.fa.*
rm data/zebraFinch/genome/Tgut_subseg_renamed.dict

rm data/BHduck/genome/MU014702.1.fa.*
rm data/BHduck/genome/MU014702.1.dict


