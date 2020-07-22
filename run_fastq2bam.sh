#!/bin/bash
#SBATCH -J sm
#SBATCH -o out
#SBATCH -e err
#SBATCH -p holy-info
#SBATCH -n 1
#SBATCH -t 4000
#SBATCH --mem=4000


#snakemake --snakefile Snakefile_fastq2bam --use-conda -p --jobs 10 --max-jobs-per-second 1 --cluster-config cluster.json --cluster "sbatch -J "sm_script" -p {cluster.p} -t {cluster.t} -n {cluster.n} -N {cluster.N} --mem={cluster.mem} "
snakemake --snakefile Snakefile_fastq2bam --profile ./profiles/slurm
