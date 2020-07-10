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

def createListsGetIndices(listDir, maxIntervalLen, maxBpPerList, maxIntervalsPerList, minNmer, ref):
    # This function ruins picard's ScatterIntervalsByNs using many MAX_TO_MERGE values to 
    # find the optimal way to balance the interval lengths (none too long) with the number
    # of intervals (not too many)

    # only create list files if they don't already exist; updating file creation dates can restart workflow!
    if len(glob.glob(listDir + "*.list")) == 0:
        command = "picard ScatterIntervalsByNs " 
        command = command + f"REFERENCE={ref} " 
        command = command + "OUTPUT=output.interval_list " 
        command = command + f"MAX_TO_MERGE={minNmer} 2> err\n"
        os.system(command)

        # raw data
        picardIntervals = defaultdict(list) # this contains the raw data to iterate back through later
        ACGTmers = defaultdict(list) 
        scaffLens = defaultdict(int)
        # info about Nmers
        NmerLens = defaultdict(int) # this just keeps track of all the unique types of Nmer lengths
        maxIntervalLenByNmer = defaultdict(int) # this stores the largest observed interval per Nmer
        NsHist = defaultdict(list) # this can be used to contruct a histogram of all the different Nmers
        # new intervals based on splitting by Nmer
        newIntervals = defaultdict(lambda: defaultdict(list))
        totalACGTmerLength = 0


        # read in picard output
        f = open("output.interval_list", 'r')
        for line in f:
            line = line.strip()
            line = line.split()
            # header of picard interval file
            if line[0].startswith("@SQ"):
                scaff = line[1].replace("SN:", "")
                length = int(line[2].replace("LN:", ""))
                scaffLens[scaff] = length 
            # store intervals, all types of Nmers in genome
            if not line[0].startswith("@"):
                scaff = line[0]
                start = int(line[1])
                end = int(line[2])
                length = end - start + 1
                merType = line[4]
                picardIntervals[scaff].append( (start, end, merType) )
                if merType == "Nmer":
                    NsHist[scaff].append(length) 
                    NmerLens[length] = 1
                elif merType == "ACGTmer":
                    ACGTmers[scaff].append( (start, end) )
                    totalACGTmerLength += length
        os.system("rm output.interval_list")
        # go through Nmers, starting with shortest ones, and construct different interval lists
        for Nmer in sorted(NmerLens):
            # go through scaffolds, from longest to shortest
            for item in sorted(scaffLens.items(), key=lambda x: x[1], reverse=True):
                scaff = item[0]
                scaffLen = item[1]
                # go through each interval on the scaffold
                currInterval = [0, 0] # use 0 to initialize, since picard does not use 0-based coords
                for i in range(len(picardIntervals[scaff])):
                    start = picardIntervals[scaff][i][0]
                    end = picardIntervals[scaff][i][1]
                    merType = picardIntervals[scaff][i][2]
                    length = end - start + 1
                    if merType == "ACGTmer":
                        # assign left coord of currInterval if it hasn't been initialized or was cleared
                        # otherwise, don't touch it!
                        if currInterval[0] == 0:
                            currInterval[0] = start
                        currInterval[1] = end
                        continue
                    if merType == "Nmer":
                        # if the Nmer is smaller than the current length we're using to split
                        # go onto the next ACGTmer and extend the interval
                        if length < Nmer:
                            continue
                        else:
                            # store currInterval bc we encounterd an Nmer of sufficient length
                            newIntervals[Nmer][scaff].append( currInterval )
                            currInterval = [0, 0]
                if currInterval[0] != 0:
                    # if the interval was extended all the way to the end of the scaffold, store it here
                    newIntervals[Nmer][scaff].append( currInterval )
        
        # for each Nmer splitting, see what the maximum interval length is
        for Nmer in sorted(NmerLens):
            totalACGTmerLength_overlaps = 0
            Nmer_maxIntervalLen = 0
            for item in sorted(scaffLens.items(), key=lambda x: x[1], reverse=True):
                scaff = item[0]
                scaffLen = item[1]
                #print(scaff, scaffLen, len(newIntervals[Nmer][scaff]))
                for itv in newIntervals[Nmer][scaff]:
                    itvLen = itv[1] - itv[0] + 1
                    if itvLen > Nmer_maxIntervalLen:
                        Nmer_maxIntervalLen = itvLen

                # below is a sanity check, to see if using this Nmer to split intervals accounts for all observed ACGTmers
                # go through new intervals, ask if they cover each of the previous ACGTmers in original picard file
                for i in ACGTmers[scaff]:
                    length = i[1] - i[0] # you usually add +1 to intervals, but this doesn't happen within the overlaps function
                    for j in newIntervals[Nmer][scaff]:
                        o = overlaps(i,j)
                        if o > 0:
                            totalACGTmerLength_overlaps += o+1 # add 1 to compare to lengths taken above

            maxIntervalLenByNmer[Nmer] = Nmer_maxIntervalLen 
            if totalACGTmerLength_overlaps < totalACGTmerLength:
                print(Nmer, totalACGTmerLength, totalACGTmerLength_overlaps)
                print("Error in interval creation: not all ACGTmers were accounted for")
                sys.exit(1)

        lessThan = []
        for Nmer in maxIntervalLenByNmer:
            if maxIntervalLenByNmer[Nmer] < maxIntervalLen:
                lessThan.append(Nmer)
                    
        if lessThan:
            print(f"best Nmer is {lessThan[-1]}")
        else:
            print(f"could not find suitable Nmer to split up genome, specify larger maximum scaffold length than {maxIntervalLen}?")
            print("Here are the maximum interval lengths we observed, for each Nmer used")
            print("Nmer MaxObservedInterval")
            for Nmer in maxIntervalLenByNmer:
                print(Nmer, maxIntervalLenByNmer[Nmer])
            sys.exit(1)

        # take optimal interval_list and use to generate GATK list files
        optimalNmer = lessThan[-1]
        runningSumBp = 0
        runningSum_intervals = 0
        current_intervals = []
        listFile_index = 0
        for scaff in newIntervals[optimalNmer]:
            for intv in newIntervals[optimalNmer][scaff]:
                start = intv[0]
                stop = intv[1]
                intervalLen = (stop - start + 1)
                if len(current_intervals) == 0:
                    # if you decide lower maxBpPerList below maxIntervalLength, it's possible tohave an uninitialized
                    # current_intervals that doesn't make it to subsequent else statement
                    current_intervals = [ (scaff, start, stop) ]
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

    # get list file indices
    LISTS = glob.glob(listDir + "*.list")	
    for i in range(len(LISTS)):
        LISTS[i] = os.path.basename(LISTS[i])
        LISTS[i] = re.search('\d+', LISTS[i]).group() # get numerical index of list
    LISTS=sorted(LISTS)
    print(LISTS)
    return(LISTS)

def overlaps(a, b):     
    """     
    Return the amount of overlap, in bp     
    between a and b.     
    If >0, the number of bp of overlap     
    If 0,  they are book-ended.     
    If <0, the distance in bp between them     
    """     
    return min(a[1], b[1]) - max(a[0], b[0])

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
                print(s, gvcfDir + s+"_L"+l+".raw.g.vcf.gz", sep="\t", file=f)
            f.close()

