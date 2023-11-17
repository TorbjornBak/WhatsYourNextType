
def Aligning(aligner, assembly, fastQFileIn):
    from progress.bar import Bar
    from Bio import Align, AlignIO,SeqIO
    
    AlignDict = dict()
    seqNr = 0
    recordSize = 0
    
    for record in SeqIO.parse(fastQFileIn,"fastq"):
        seqNr += 1
        if recordSize < len(record.seq) and len(record.seq) < 4000:
            firstRecord = record
            recordSize = len(firstRecord.seq)
    
    bar = Bar('Aligning', max = seqNr, fill = "|")
    cutoff = 3000

    Alignments_tsv = open("Alignments2_tsv", "w")
    firstTime = True
    ReversedFastq = open("HLAC_bin.fastq","w")
    for record in SeqIO.parse(fastQFileIn,"fastq"):
        if len(record.seq) > cutoff:
            for assemblies in SeqIO.parse(assembly,"fasta"):
                #Use assembly
                firstRecord = assemblies
                
                length = min(len(firstRecord.seq), len(record.seq))
                forwardAlign =aligner.align(firstRecord.seq[:length],record.seq[2:length-3])
                revAlign = aligner.align(firstRecord.seq[:length],record.reverse_complement()[2:length-3])
                if firstTime == True:
                    print("Reference", "\t", firstRecord.seq, file = Alignments_tsv)
                    firstTime = False
                if forwardAlign.score > revAlign.score:
                    alignment = forwardAlign[0]
                    break
                else:
                    alignment = revAlign[0]
                    break
                    
            print("@"+record.id, "\t", alignment[1,:], file = Alignments_tsv)
            print()
            AlignDict["@" + record.id] = alignment[1,:]
        bar.next()
    bar.finish()
    return AlignDict

def fromTSV(fileName):
    TSV_file = open(fileName)
    AlignDict = dict()
    for line in TSV_file:
        id = line.split("\t")[0].strip()
        seq = line.split("\t")[1].strip()
        AlignDict[id] = seq
    return AlignDict

def Frequency(AlignDict):
    FreqDict = dict()
    for key,value in AlignDict.items():
        for i in range(len(value)):
            if str(i) in FreqDict:
                if value[i] in FreqDict[str(i)]:
                    FreqDict[str(i)][value[i]] += 1
                else:
                    FreqDict[str(i)][value[i]] = 1
            else:
                FreqDict[str(i)] = dict()
    return FreqDict


def plotit(FreqDict):
    import matplotlib.pyplot as plt
    YList = []
    AList = []
    GList = []
    CList = []
    TList = []
    NList = []
    DashList = []
    #print(FreqDict)
    for key,value in FreqDict.items():
        YList.append(key)
        if "C" in value.keys():
            CList.append(value["C"])
        else:
            CList.append(0)
            
        if "G" in value.keys():
            GList.append(value["G"])
        else:
            GList.append(0)
            
        if "A" in value.keys():
            AList.append(value["A"])
        else:
            AList.append(0)
            
        if "T" in value.keys():
            TList.append(value["T"])
        else:
            TList.append(0)
            
        if "N" in value.keys():
            NList.append(value["N"])
        else:
            NList.append(0)
            
        if "-" in value.keys():
            DashList.append(value["-"])
        else:
            DashList.append(0)
            
    plt.figure(figsize=(50,10))
    ListList = [AList,GList,CList,TList, NList]
    BaseList = ["A","G","C","T","N"]
    colours = ["green","blue","yellow","red","black"]
    for i in range(len(ListList)):
        if len(YList) > len(ListList[i]):
            for p in range(len(YList)-len(ListList[i])):
             ListList[i].append(0)
        plt.subplot(231 + i)
        plt.plot(YList,ListList[i][:len(YList)], c=colours[i])
        plt.grid(True)
        plt.title(BaseList[i], fontsize = 30)
    plt.savefig("Densities.png")
    
    nMean = sum(NList)/sum([sum(CList),sum(GList),sum(NList),sum(AList),sum(TList)])
    return nMean

def main(assembly,fastQFileIn):
    import numpy as np
    import pandas as pd

    from Bio import Align, AlignIO,SeqIO
    
    import sys
    
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"

    aligner.target_left_extend_gap_score = -100
    aligner.target_left_open_gap_score =  -100
    aligner.target_right_open_gap_score =  -100
    aligner.target_right_extend_gap_score = -100
    aligner.target_internal_open_gap_score = -100
    aligner.target_internal_extend_gap_score =-20


    aligner.match_score = 3
    aligner.mismatch_score = -2

    aligner.query_right_extend_gap_score = -0.5
    aligner.query_right_open_gap_score = -1
    aligner.query_left_open_gap_score = -1
    aligner.query_left_extend_gap_score = -0.5
    aligner.query_internal_open_gap_score = -1
    aligner.query_internal_extend_gap_score = -0.5
    
    AlignDict = Aligning(aligner,assembly, fastQFileIn)
    FreqDict = Frequency(AlignDict)
    
    return FreqDict


