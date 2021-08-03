#!/usr/bin/python -tt
import glob
import os
import sys
import re
from collections import defaultdict
from urllib.request import urlopen

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

def printBamSumStats(FractionReadsPassFilter, NumFilteredReads, PercentDuplicates, PercentHQreads, PercentHQbases, SeqDepths, CoveredBases, validateSams, outputDir, wildcards):
    out_path = os.path.join(outputDir, wildcards.Organism, wildcards.refGenome, "bam_sumstats.txt")
    o = open(out_path, 'w')
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

def getFastqSampleNames(fastqDir, fastq_suffix1):
    SAMPLES = glob.glob(fastqDir + "*" + fastq_suffix1)	
    for i in range(len(SAMPLES)):
        SAMPLES[i] = os.path.basename(SAMPLES[i])
        SAMPLES[i] = SAMPLES[i].replace(fastq_suffix1, "")
    #print(SAMPLES)
    return(SAMPLES)
    
def getBamSampleNames(BamDir, bam_suffix):
    SAMPLES = glob.glob(BamDir + "*" + bam_suffix)	
    for i in range(len(SAMPLES)):
        SAMPLES[i] = os.path.basename(SAMPLES[i])
        SAMPLES[i] = SAMPLES[i].replace(bam_suffix, "")
    return(SAMPLES)

def createSeqDictGetScaffOrder(dict_file):
    #print(f"{ref=}, {refBaseName=}")
    #seqDict = refBaseName + ".dict"
    #index = ref + ".fai"
    seqDictScaffs = []
    f = open(dict_file, 'r')
    for line in f:
        if line.startswith("@SQ"):
            line = line.split("\t")
            scaff = line[1].replace("SN:", "")
            seqDictScaffs.append(scaff)
    return(seqDictScaffs)

def createListsGetIndices(maxIntervalLen, maxBpPerList, maxIntervalsPerList, minNmer, outputDir, intDir, wildcards, dict_file, intervals_file):

    # Get the order in which scaffolds are listed in genome .dict file,  need to correspond to order in list files! So use .dict order for printing interval files
    seqDictScaffs = createSeqDictGetScaffOrder(dict_file)

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
    
    #get actual basename. refBaseName is defined earlier and its usage requires the full path
    #refFileName = refBaseName.rsplit("/")[-1]

    # read in picard output
    f = open(intervals_file, 'r')
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

    # if picard could not find any Nmers in your assembly (possible if only looking at small region) give NmerLens a key of 1 so loop below executes
    if len(NmerLens.keys()) == 0:
        NmerLens[0] = 1
    
    # go through Nmers, starting with shortest ones, and construct different interval lists
    # but skip Nmers within 50bp of one another as they probably give similar results.
    prevNmer = 0
    for Nmer in sorted(NmerLens):
        # need to initialize with prevNmer == 1 in case there are no Nmers found, such that the only one is of length one (see above code)
        if prevNmer == 0 or Nmer >= (prevNmer + 50):
            prevNmer = Nmer
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
    for Nmer in sorted(newIntervals.keys()):
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
                
    outfile = os.path.join(outputDir, wildcards.Organism, wildcards.refGenome, intDir, wildcards.refGenome + "_interval_algo.out")
    out = open(outfile,'w')
    print("Here are the actual maximum interval lengths we observed, for each minimum Nmer size used to split up the genome.", file=out)
    print("To proceed, specify a maxIntervalLen that is equal to or slightly greater than the ones observed here.", file=out)
    print("Nmer MaxObservedInterval NumIntervals", file=out)
    for Nmer in maxIntervalLenByNmer:
        numInt = 0
        for scaff in newIntervals[Nmer]:
            numInt += len(newIntervals[Nmer][scaff])
        print(Nmer, maxIntervalLenByNmer[Nmer], numInt, file=out)
    out.close()

    if lessThan:
        print(f"best Nmer is {lessThan[-1]}")
    else:
        print(f"We could not find suitable Nmer to split up genome given the specified maxIntervalLen (which is currently too small).")
        print("Here are the actual maximum interval lengths we observed, for each minimum Nmer size used to split up the genome.")
        print("To proceed, specify a maxIntervalLen that is equal to or slightly greater than the ones observed here.")
        print("Nmer MaxObservedInterval NumIntervals")
        for Nmer in maxIntervalLenByNmer:
            numInt = 0
            for scaff in newIntervals[Nmer]:
                numInt += len(newIntervals[Nmer][scaff])
            print(Nmer, maxIntervalLenByNmer[Nmer], numInt)
        sys.exit("\n\n\n*****PLEASE SEE out FILE FOR A MESSAGE ON HOW TO PROCEED! :)*****\n\n\n")

    # take optimal interval_list and use to generate GATK list files
    optimalNmer = lessThan[-1]
    runningSumBp = 0
    runningSum_intervals = 0
    current_intervals = []
    listFile_index = 0
    # go through scaffolds in same order as they appear in the sequence dictionary, printing lists as we go
    outfile = os.path.join(outputDir, wildcards.Organism, wildcards.refGenome, intDir, wildcards.refGenome + "_intervals_fb.bed")
    bed = open(outfile, 'w')
    for scaff in seqDictScaffs:
        for intv in newIntervals[optimalNmer][scaff]:
            start = intv[0]
            stop = intv[1]
            intervalLen = (stop - start + 1)
            print(scaff, start, stop, sep="\t", file=bed)
            # does adding the next interval put us over the maximum values for our two thresholds?
            if (runningSumBp + intervalLen) >= maxBpPerList or (runningSum_intervals + 1) >= maxIntervalsPerList:
                if len(current_intervals) == 0:
                    # if you decide lower maxBpPerList below maxIntervalLength, it's possible tohave an uninitialized
                    # current_intervals that doesn't make it to subsequent else statement
                    current_intervals = [ (scaff, start, stop) ]
                # flush out current_intervals into a list file
                printIntervalsToListFile(intDir, listFile_index, current_intervals, wildcards, outputDir)
                #out = open(f"{intDir}list{listFile_index}.list", 'w')
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
        printIntervalsToListFile(intDir, listFile_index, current_intervals, wildcards, outputDir)
        #out = open(f"{intDir}list{listFile_index}.list", 'w')
        #for i in current_intervals:
        #    print(f"{i[0]}:{i[1]}-{i[2]}", file=out)
        #out.close()
    bed.close()

    # get list file indices
    gatk_list_dir = os.path.join(outputDir, wildcards.Organism, wildcards.refGenome, intDir)
    LISTS = glob.glob(gatk_list_dir + "/*.list")
    print(LISTS)	
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

def printIntervalsToListFile(intDir, listFile_index, current_intervals, wildcards, outputDir):
    gatk_list_dir = os.path.join(outputDir, wildcards.Organism, wildcards.refGenome, intDir)
    
    if not os.path.isdir(gatk_list_dir):
        os.mkdir(gatk_list_dir)
    
    outfile = os.path.join(gatk_list_dir, f"list{listFile_index}.list")
    out = open(outfile, 'w')
    
    for i in current_intervals:
        print(f"{i[0]}:{i[1]}-{i[2]}", file=out)
    out.close()


def loadIntervalsForFB(f):
    intervals_fb = []
    fh = open(f, 'r')
    for line in fh:
        line = line.split()
        interval = line[0] + ":" + line[1] + "-" + line[2]
        intervals_fb.append(interval)

    return(intervals_fb)
    
def getVcfs_gatk(LISTS, vcfDir):
    vcfs = []
    for i in range(len(LISTS)):
       vcfs.append(f"-I {vcfDir}L{i}_filter.vcf")
    out = " ".join(vcfs)
    return(out)

def getListIndices(intDir):
    LISTS = glob.glob(intDir + "gatkLists/*.list")	
    for i in range(len(LISTS)):
        LISTS[i] = os.path.basename(LISTS[i])
        LISTS[i] = re.search('\d+', LISTS[i]).group() # get numerical index of list
    LISTS=sorted(LISTS)
    return(LISTS)


