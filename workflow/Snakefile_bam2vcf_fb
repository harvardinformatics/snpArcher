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

# this is where Snakemake output will go, specify with baseDir in config.yml
fbDir = config["fbDir"]
bamDir = config["bamsForFB"]
intervalDir = fbDir + config["intervalDir"]  
vcfDir_fb = fbDir + config["vcfDir_fb"]
intDir = config["intDir"]
maxDP_fb = config["maxDP_fb"]

refBaseName = helperFun.getRefBaseName(config["ref"])

# grab all samples for R1 to get list of names, no need to look at R2 which should have identical names
#SAMPLES = ["ERR1013163"]
SAMPLES = helperFun.getBamSampleNames(bamDir, bam_suffix)


if not os.path.isdir(fbDir):
    os.system("mkdir " + fbDir)
intervals_fb = helperFun.loadIntervalsForFB(intDir + "intervals_fb.bed")


def myfunc(vcfDir_fb, intervals_fb):
    gatherCommand = ""
    for j in intervals_fb:
        gatherCommand = gatherCommand + f"-I {vcfDir_fb}{j}.vcf "
    return(gatherCommand)

### workflow ###

rule all:
    input:
        missing = config["fbDir"] + "missing_data_per_ind.txt",
        SNPsPerInt = config["fbDir"] + "SNP_per_interval.txt"

include: "rules/bam2vcf_fb.smk"

