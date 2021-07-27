import os
from collections import defaultdict
import sys
import yaml
sys.path.append(os.getcwd())
import helperFun

configfile: "config.yaml"
res_config = yaml.load(open("resources.yaml"))

# rename variables from config file for downstream simplicity
#ref = config["ref"]
bam_suffix = config["bam_suffix"]

# this is where Snakemake output will go, specify with 'gatkDir' in config.yml
bamDir = config["bamsForGatk"] 
intDir = config["intDir"]
maxIntervalLen = int(config["maxIntervalLen"])
maxBpPerList = int(config["maxBpPerList"])
maxIntervalsPerList = int(config["maxIntervalsPerList"])
minNmer = config["minNmer"]

#refBaseName = helperFun.getRefBaseName(config["ref"])

# grab all samples for R1 to get list of names, no need to look at R2 which should have identical names
SAMPLES = helperFun.getBamSampleNames(bamDir, bam_suffix)

# create directory to store GATK lists
if not os.path.isdir(intDir + "gatkLists"):
    os.system("mkdir -p " + intDir + "gatkLists")

_, sample_dict, _ = helperFun.create_sample_dict(config["samples"])

GENOMES = {sample_dict[k]["refGenome"]: sample_dict[k]["Organism"] for k in sample_dict.keys()}
print(GENOMES)
### workflow ###

rule all:
    input:
        intervals = expand(intDir + "{genome}_intervals_fb.bed", genome=GENOMES.keys()),


rule process_ref:
    """
    This rule generates a .fai file from the reference genome, which are required for GATK to run properly. GATK also needs a .dict file, but this was previously generated.
    """
    input:
        ref = lambda wildcards: 
                    expand("data/{organism}/genome/{{genome}}.fna", organism=GENOMES[wildcards.genome])
    output: 
        fai = "data/{organism}/genome/{genome}.fna.fai",
        dictf = "data/{organism}/genome/{genome}.dict",
    conda:
        "./envs/bam2vcf.yml"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * res_config['process_ref']['mem']   
    shell:
        "samtools faidx {input.ref} --output {output.fai}\n"
        "picard CreateSequenceDictionary REFERENCE={input.ref} OUTPUT={output.dictf}\n"
        
rule create_intervals:
    input:
        fai = lambda wildcards: 
                expand("data/{organism}/genome/{{genome}}.fna.fai", organism=GENOMES[wildcards.genome]),
        dictf = lambda wildcards: 
                expand("data/{organism}/genome/{{genome}}.dict", organism=GENOMES[wildcards.genome]),
        ref = lambda wildcards: 
                    expand("data/{organism}/genome/{{genome}}.fna", organism=GENOMES[wildcards.genome]),
        intervals = intDir + "{genome}_output.interval_list"
    params:
        refBaseName = lambda wildcards: 
                expand("data/{organism}/genome/{genome}", genome=wildcards.genome, organism=GENOMES[wildcards.genome]),
    output: 
        intervals = intDir + "{genome}_intervals_fb.bed"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * res_config['create_intervals']['mem'] 
    run:
        LISTS = helperFun.createListsGetIndices(intDir, maxIntervalLen, maxBpPerList, maxIntervalsPerList, minNmer, input.ref[0], params.refBaseName[0], intDir)

rule picard_intervals:
    input:
        ref = lambda wildcards: 
                expand("data/{organism}/genome/{{genome}}.fna", organism=GENOMES[wildcards.genome])
    output:
        intervals = temp(intDir + "{genome}_output.interval_list")
    conda:
        "./envs/bam2vcf.yml"
    log:
        "log/picard_intervals/{genome}.log"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * res_config['process_ref']['mem']   
    shell:
        "picard ScatterIntervalsByNs REFERENCE={input} OUTPUT={output.intervals} MAX_TO_MERGE={minNmer} > {log}\n" 



