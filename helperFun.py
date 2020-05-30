#!/usr/bin/python -tt
import glob
import os
import re
from collections import defaultdict

def collectDedupMetrics(dedupFiles):

    PercentDuplicates = defaultdict(float)
    for fn in dedupFiles:
        sample = os.path.basename(fn)
        sample = sample.replace("_dedupMetrics.txt", "")
        f = open(fn, 'r')
        lines = f.readlines()
        for i in range(len(lines)):
            if lines[i].startswith("LIBRARY"):
                info = lines[i+1]
                info = info.split("\t")
                percDup = info[8]
                #print(percDup)
                PercentDuplicates[sample] = percDup
        f.close()
    return(PercentDuplicates)

def collectAlnSumMets(alnSumMetsFiles):

    PercentHQreads = defaultdict(float)
    PercentHQbases = defaultdict(float)

    for fn in alnSumMetsFiles:
        sample = os.path.basename(fn)
        sample = sample.replace("_AlnSumMets.txt", "")
        f = open(fn, 'r')
        for line in f:
            if line.startswith("PAIR"):
                line = line.split()
                totalReads = int(line[1])
                totalAlignedBases = int(line[7])
                HQreads = int(line[8])
                totalAlignedHQ20bases = int(line[10])
                PercentHQreads[sample] = HQreads/totalReads
                PercentHQbases[sample] = totalAlignedHQ20bases/totalAlignedBases
        f.close()
    return(PercentHQreads, PercentHQbases)

def collectValidationStatus(validateFiles):

    validateSams = defaultdict(float)

    for fn in validateFiles:
        sample = os.path.basename(fn)
        sample = sample.replace("_validate.txt", "")
        f = open(fn, 'r')
        status = "FALSE"
        for line in f:
            if "No errors found" in line:
                status = "TRUE"
        f.close()
        validateSams[sample] = status
    return(validateSams)

def collectCoverageMetrics(coverageFiles):

    SeqDepths = defaultdict(float)
    CoveredBases = defaultdict(float)

    for fn in coverageFiles:
        # these files contain coverage data by scaffold; take weighted average
        sample = os.path.basename(fn)
        sample = sample.replace("_coverage.txt", "")
        numSites = []
        covbases = 0
        depths = [] # samtools coverage function prints depth per scaffold
        f = open(fn, 'r')
        for line in f:
            if not line.startswith("#rname"):
                line = line.split()
                numSites.append( int(line[2]) - int(line[1]) + 1 )
                depths.append( float(line[6]) )
                covbases +=  float(line[4])
        f.close()
        total = sum(numSites)
        depthsMean = 0
        for i in range(len(depths)):
            depthsMean += depths[i]*numSites[i]/(total)
        SeqDepths[sample] = depthsMean
        CoveredBases[sample] = covbases
    return(SeqDepths, CoveredBases)

def printBamSumStats(PercentDuplicates, PercentHQreads, PercentHQbases, SeqDepths, CoveredBases, validateSams):

    o = open("bam_sumstats.txt", 'w')
    print("sample", "PercentDuplicates", "PercentHQ20alignedReads", "PercentHQ20bases", "MeanSeqDepth", "BasesCoveredMoreThanOnce", "validBAM", file=o, sep="\t")
    for sample in PercentDuplicates:
        print(sample,file=o, end="\t")
        print(PercentDuplicates[sample], file=o, end="\t")
        print(PercentHQreads[sample], file=o, end="\t")
        print(PercentHQbases[sample], file=o, end="\t")
        print(SeqDepths[sample], file=o, end="\t")
        print(CoveredBases[sample], file=o, end="\t")
        print(validateSams[sample], file=o)

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
        fileName = dbDir + "DB_mapfile" + l
        # only rewrite file if you haven't already; otherwise rerunning the script updates time stamp of these files
        # causing snakemake to erroneously rerun some steps due to updated input
        if not os.path.exists(fileName):
            f=open(fileName, 'w')
            for s in SAMPLES:
                print(s, gvcfDir + s+"_L"+l+".raw.g.vcf", sep="\t", file=f)
            f.close()

