import numpy as np
import pandas as pd
import sys


from Bio import pairwise2
from Bio import SeqIO




assembly = sys.argv[1]
fileName = sys.argv[2]
fastQFileIn = sys.argv[3]
path = sys.argv[4]

DistDict1 = dict()
ClusterDict = dict()
#Check the file to see if we have 1 or more contigs
twoFiles = 0
for prerecord in SeqIO.parse(assembly,"fasta"):
    twoFiles += 1

if twoFiles == 1:
    for record in SeqIO.parse(fastQFileIn,"fastq"):
        for assemblies in SeqIO.parse(assembly,"fasta"):
            forwardAlign = pairwise2.align.globalxx(record.seq,assemblies.seq, score_only = True)
            revAlign = pairwise2.align.globalxx(record.reverse_complement(),assemblies.seq, score_only = True)
            bestAlignment = max(forwardAlign,revAlign)
            DistDict1["@"+record.id] = bestAlignment
    
    avgDist = np.mean(np.array(list(DistDict1.values())))
    Cluster1 =[key for key,value in DistDict1.items() if (value) > avgDist]
    Cluster2 =[key for key,value in DistDict1.items() if (value) <= avgDist]
    FastQFile = open(fastQFileIn,"r")
    OutputFile1 = open(path + "/Clusters_{}/".format(fileName) + fileName+".cluster1.fastq","w")
    OutputFile2 = open(path + "/Clusters_{}/".format(fileName)  + fileName+".cluster2.fastq","w")

    FastaFlag1 = False
    FastaFlag2 = False
    iterator = 0

    for record in SeqIO.parse(fastQFileIn,"fastq"):
        if "@"+record.id in Cluster1:
            print(record.format("fastq"), file = OutputFile1, end="")
        if "@"+record.id in Cluster2:
            print(record.format("fastq"), file = OutputFile2, end="")
            

else:
    for assemblies in SeqIO.parse(assembly,"fasta"):
        recordnumber = 0
        bestAlignment = 0
        for record in SeqIO.parse(fastQFileIn,"fastq"):        
            recordnumber += 1
            forwardAlign = pairwise2.align.globalxx(record.seq,assemblies.seq, score_only = True)
            revAlign = pairwise2.align.globalxx(record.reverse_complement(),assemblies.seq, score_only = True)
            if max(forwardAlign, revAlign) > bestAlignment:
                bestAlignment = max(forwardAlign,revAlign)
                DistDict1["@"+record.id] = [bestAlignment,recordnumber]
    for i in range(twoFiles):
        ClusterDict["cluster"+str(i)] = [key for key,value in DistDict1.items() if value[1] == i]
    for Clusters in ClusterDict.keys():
        OutputFile = open(path + "/Clusters_{}/".format(fileName) + fileName+".{}.fastq".format(Clusters),"w")
        for record in SeqIO.parse(fastQFileIn,"fastq"):
            if "@" + record.id in ClusterDict[Clusters]:
                print(record.format("fastq"), file = OutputFile, end="")


    
    





