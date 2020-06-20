#!/usr/bin/python -tt
import glob
import os
import sys
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

def collectFastpOutput(fastpFiles):

    FractionReadsPassFilter = defaultdict(float)
    NumFilteredReads = defaultdict(int)

    for fn in fastpFiles:
        sample = os.path.basename(fn)
        sample = sample.replace("_fastp.out", "")
        unfiltered = 0
        filtered = 0
        f = open(fn, 'r')
        for line in f:
            if "before filtering" in line:
                line = next(f)
                line = line.split()
                unfiltered += int(line[2])
            if "Filtering result" in line:
                line = next(f)
                line = line.split()
                filtered = int(line[3])
        f.close()
        FractionReadsPassFilter[sample] = float(filtered/unfiltered)
        NumFilteredReads[sample] = filtered
        
    return(FractionReadsPassFilter, NumFilteredReads)

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

def printBamSumStats(FractionReadsPassFilter, NumFilteredReads, PercentDuplicates, PercentHQreads, PercentHQbases, SeqDepths, CoveredBases, validateSams):

    o = open("bam_sumstats.txt", 'w')
    print("sample", "FractionReadsPassFilter", "NumFilteredReads", "PercentDuplicates", "PercentHQ20alignedReads", "PercentHQ20bases", "MeanSeqDepth", "BasesCoveredMoreThanOnce", "validBAM", file=o, sep="\t")
    for sample in PercentDuplicates:
        print(sample,file=o, end="\t")
        print(FractionReadsPassFilter[sample], file=o, end="\t")
        print(NumFilteredReads[sample], file=o, end="\t")
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
    elif ".fna" in ref:
        refBaseName = ref.replace(".fna", "")
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

def createListsGetIndices(listDir, maxIntervalLen, maxBpPerList, maxIntervalsPerList, ref):
    # This function ruins picard's ScatterIntervalsByNs using many MAX_TO_MERGE values to 
    # find the optimal way to balance the interval lengths (none too long) with the number
    # of intervals (not too many)

    # only create list files if they don't already exist; updating file creation dates can restart workflow!
    if len(glob.glob(listDir + "*.list")) == 0:
        NsPassed = []
        maxLengthObserved = 0
        # run picard's ScatterIntervalsByNs using many MAX_TO_MERGE values
        for Ns in range(10000, 150000, 10000):
            command = "picard ScatterIntervalsByNs " 
            command = command + f"REFERENCE={ref} " 
            command = command + f"OUTPUT_TYPE=ACGT OUTPUT=Ns{Ns}.interval_list "
            command = command + f"MAX_TO_MERGE={Ns} 2> out\n"
            os.system(command)
        # analyze picard output, recording the maximum interval length for each scattering
        for Ns in range(10000, 150000, 10000):
            f = open(f"Ns{Ns}.interval_list", 'r')
            for line in f:
                line = line.split()
                if not line[0].startswith("@"):
                    length = int(line[2]) - int(line[1])
                    if length > maxLengthObserved:
                        maxLengthObserved = length
            if maxLengthObserved <= maxIntervalLen:
                NsPassed.append(Ns)
            f.close()
        if NsPassed:
            print(NsPassed[-1])
        else:
            sys.exit("error in interval files creation")

        # take optimal interval_list and use to generate GATK list files
        NsForSplitting = NsPassed[-1]
        runningSumBp = 0
        runningSum_intervals = 0
        current_intervals = []
        listFile_index = 0
        f = open(f"Ns{NsForSplitting}.interval_list", 'r')
        for line in f:
            line = line.split()
            # skip file header, starts with @
            if not line[0].startswith("@"):
                scaff = line[0]
                start = int(line[1])
                stop = int(line[2])
                intervalLen = (stop - start)
                # does adding the next interval put us over the maximum values for our two thresholds?
                if (runningSumBp + intervalLen) >= maxBpPerList or (runningSum_intervals + 1) >= maxIntervalsPerList:
                    # flush out current_intervals into a list file
                    printIntervalsToListFile(listDir, listFile_index, current_intervals)
                    #out = open(f"{listDir}list{listFile_index}.list", 'w')
                    #for i in current_intervals:
                    #    print(f"{i[0]}:{i[1]}-{i[2]}", file=out)
                    #out.close()
                    # re-initialize data for next list file
                    current_intervals = [ (scaff, start, stop) ]
                    runningSumBp = intervalLen
                    runningSum_intervals = 1 
                    listFile_index += 1
                else:
                    current_intervals.append( (scaff, start, stop) )
                    runningSum_intervals += 1
                    runningSumBp += intervalLen
        # if part-way through the loop above ran out of intervals to print, print whatever remains
        if current_intervals:
            printIntervalsToListFile(listDir, listFile_index, current_intervals)
            #out = open(f"{listDir}list{listFile_index}.list", 'w')
            #for i in current_intervals:
            #    print(f"{i[0]}:{i[1]}-{i[2]}", file=out)
            #out.close()
        os.system("rm Ns*.interval_list")

    # get list file indices
    LISTS = glob.glob(listDir + "*.list")	
    for i in range(len(LISTS)):
        LISTS[i] = os.path.basename(LISTS[i])
        LISTS[i] = re.search('\d+', LISTS[i]).group() # get numerical index of list
    LISTS=sorted(LISTS)
    print(LISTS)
    return(LISTS)

def printIntervalsToListFile(listDir, listFile_index, current_intervals):
    out = open(f"{listDir}list{listFile_index}.list", 'w')
    for i in current_intervals:
        print(f"{i[0]}:{i[1]}-{i[2]}", file=out)
    out.close()

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

