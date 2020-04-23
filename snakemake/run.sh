#!/bin/bash
#SBATCH -J sm
#SBATCH -o out2
#SBATCH -e err2
#SBATCH -p holy-info
#SBATCH -n 1
#SBATCH -t 4000
#SBATCH --mem=4000


#snakemake --dryrun --jobs 3 --cluster-config cluster.json --cluster "sbatch -J "test" -p {cluster.p} -t {cluster.t} -n {cluster.n} -N {cluster.N} --mem={cluster.mem} "
snakemake -p --jobs 300 --cluster-config cluster.json --cluster "sbatch -J "sm_script" -p {cluster.p} -t {cluster.t} -n {cluster.n} -N {cluster.N} --mem={cluster.mem} "

