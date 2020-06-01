#!/usr/bin/python -tt

import re
import sys
import os
from Bio import SeqIO
from collections import defaultdict

def getLengths(fileName):
    scaffLengths = []
    totalLen = 0
    
    f = open(fileName)
    for record in SeqIO.parse(fileName, "fasta"):
        scaff = record.id
        length = len(record.seq)
        scaffLengths.append( (scaff, length) )
        totalLen += length
    return(sorted(scaffLengths, key = lambda x: x[1], reverse=True), totalLen)

def createListsDict(scaffLengths, segmentLens, maxScaffsPerBin):
    listsDict = defaultdict(list)
    currCumulativeLen = 0
    currList = []
    currBin = 0
    for i in range(len(scaffLengths)): 
        scaff = scaffLengths[i][0]
        length = scaffLengths[i][1]
        currCumulativeLen += length
        currList.append(scaff)

        # if you've gone over sequence length for a bin or the max number of scaffolds per bin
        # store current contents and reset
        if (int((currCumulativeLen)/segmentLens) >= 1) or (len(currList) >= maxScaffsPerBin):
            listsDict[currBin] = currList
            currList = []
            currCumulativeLen = 0
            currBin += 1
        else:
            next
    if len(currList) > 0: 
        # if you reached the end of the fasta file and did not meet the contidions above
        # store the currList in the currBin
        listsDict[currBin] = currList
    return(listsDict)

def main():

    fileName = sys.argv[1]
    numGenomeSegments = 20 # this is the number of bins to divide the scaffolds into, but the last bin is treated differently since it contains too many, goes too slow
    maxScaffsPerBin = 200 # for all the scaffolds remaining in the last bin, this divides them up by this number into additional bins 

    # get the length of each scaffold (scaffLengths) and the total length of the genome (totalLen)
    scaffLengths, totalLen = getLengths(fileName)
    print( "Number of scaffolds:", len(scaffLengths))
    print( "Total genome length", totalLen)

    # calculate the minimum number of bp contained in each list (that may consist of one or more scaffolds)
    segmentLens = int(totalLen/numGenomeSegments) 

    listsDict = createListsDict(scaffLengths, segmentLens, maxScaffsPerBin)

    for i in sorted(listsDict.keys()):
        outFile = open('list%s.list' % i, 'w')
        for scaff in listsDict[i]:
            print(scaff, file=outFile)

if __name__ == '__main__':
  main()
