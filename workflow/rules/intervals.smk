import os
from collections import defaultdict
import sys
import yaml
sys.path.append(os.getcwd())
import helperFun



configfile: "config.yaml"
res_config = yaml.load(open("resources.yaml"))

# rename variables from config file for downstream simplicity
ref = config["ref"]
bam_suffix = config["bam_suffix"]

# this is where Snakemake output will go, specify with 'gatkDir' in config.yml
bamDir = config["bamsForGatk"] 
intDir = config["intDir"]

# create directory to store GATK lists
if not os.path.isdir(intDir + "gatkLists"):
    os.system("mkdir -p " + intDir + "gatkLists")

_, sample_dict, _ = helperFun.create_sample_dict("samples.csv")

GENOMES = {sample_dict[k]["refGenome"]: sample_dict[k]["Organism"] for k in sample_dict.keys()}
print(GENOMES)
print(GENOMES.keys())

### workflow ###

rule all:
    input: expand("data/Escherichia_coli/genome/{genome}{ext}/", genome=GENOMES.keys(), ext=[".fai", ".dict"])

rule process_ref:
    """
    This rule generates a .fai file from the reference genome, which are required for GATK to run properly. GATK also needs a .dict file, but this was previously generated.
    """
    input:
        ref = lambda wildcards: 
            expand("data/{organism}/genome/{{genome}}.fna", organism=GENOMES[wildcards.genome])
    output: 
        fai = "data/{organism}/genome/{genome}.fai",
        dictf = "data/{organism}/genome/{genome}.dict",
    conda:
        "./envs/bam2vcf.yml"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * res_config['process_ref']['mem']   
    shell:
        "samtools faidx {input.ref} --output {output.fai}\n"
        "picard CreateSequenceDictionary REFERENCE={input.ref} OUTPUT={output.dictf}\n"