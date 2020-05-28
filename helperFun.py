#!/usr/bin/python -tt
import glob
import os
import re

def getRefBaseName(ref):
    if ".fasta" in ref:
        refBaseName = ref.replace(".fasta", "")
    elif ".fa" in ref:
        refBaseName = ref.replace(".fa", "")
    else:
        sys.exit("Your reference genome did not end in \".fasta\" or \".fa\". Please fix this before continuing.")
    return(refBaseName)

def getSampleNames(fastqDir, fastq_suffix1):
    SAMPLES = glob.glob(fastqDir + "*" + fastq_suffix1)	
    for i in range(len(SAMPLES)):
        SAMPLES[i] = os.path.basename(SAMPLES[i])
        SAMPLES[i] = SAMPLES[i].replace(fastq_suffix1, "")
    print(SAMPLES)
    return(SAMPLES)

def getListIndices(listDir):
    LISTS = glob.glob(listDir + "*.list")	
    for i in range(len(LISTS)):
        LISTS[i] = os.path.basename(LISTS[i])
        LISTS[i] = re.search('\d+', LISTS[i]).group() # get numerical index of list
    LISTS=sorted(LISTS)
    print(LISTS)
    return(LISTS)

def makeMapFilesForGenomicsDBImport(SAMPLES, LISTS, dbDir, gvcfDir):
    for l in LISTS:
        f=open(dbDir + "DB_mapfile" + l, 'w')
        for s in SAMPLES:
            print(s, gvcfDir + s+"_L"+l+".raw.g.vcf", sep="\t", file=f)
        f.close()

