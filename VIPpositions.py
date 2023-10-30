import sys

def VIPPos(FreqDict):
    VIPDict = dict()
    for key1 in FreqDict.keys():
        if len(FreqDict[key1]) < 1:
            break
        SortedFreq = sorted(FreqDict[key1].items(), key=lambda x:x[1],reverse=True)
        VIPDict[key1] = SortedFreq[0:2]
        
    return VIPDict

def NReplace(VIPDict, fastqfile, tsvFile ):
    from Bio import SeqIO
    from Clustering import Frequency
    
    qualityDict = []
    
    for record in SeqIO.parse(fastqfile,"fastq"):
        if len(record.seq) > 3000:
            qualityDict.append(("@"+record.id, record.letter_annotations["phred_quality"]))
    
    
    alignFile = open(tsvFile, "r")
    outputFastq = open("N"+fastqfile[-14:], "w")
    NDict = dict()
    recordIterator = 0
    
    for line in alignFile:
        recordID = line.split("\t")[0].strip()
        seq = line.split("\t")[1]
        i = 0
        newSeq= []
        NDict[recordID] = ""
        if recordID != "Reference":
            print(recordID, file = outputFastq)
            QualityLine = ""
            QualityIterator = 0
            print(recordID)
            print(len(seq))
            print(qualityDict[recordIterator][0])
            print(len(qualityDict[recordIterator][1]))
            
            for char in seq[:len(VIPDict.keys())-6]:
                if QualityIterator > len(seq):
                    break
                if len(VIPDict[str(i)]) > 1:
                    if char == "-":
                        print("N",file = outputFastq, end = "")
                        QualityLine += chr(33)
                        
                    elif char == VIPDict[str(i)][0][0]:
                        print(char,file = outputFastq, end = "")
                        QualityLine += chr(qualityDict[recordIterator][1][QualityIterator]+33)
                        QualityIterator += 1
                        
                    elif char == VIPDict[str(i)][1][0]:
                        print(char,file = outputFastq, end = "")
                        QualityLine += chr(qualityDict[recordIterator][1][QualityIterator]+33)
                        QualityIterator += 1
                        
                    else:
                        print("N",file = outputFastq, end = "")
                        QualityLine += chr(33)
                        
                else:
                    if char == "-":
                        print("N",file = outputFastq, end = "")
                        QualityLine += chr(33)
                        
                    elif char == VIPDict[str(i)][0][0]:
                        print(char,file = outputFastq, end = "")
                        QualityLine += chr(qualityDict[recordIterator][1][QualityIterator]+33)
                        QualityIterator += 1
                        
                    else:
                        print("N",file = outputFastq, end = "")
                        QualityLine += chr(33)
                
                 
                i+= 1
            print("\n+", file = outputFastq)
            print(QualityLine, file = outputFastq)
            recordIterator += 1
    FreqDict2 = Frequency(NDict)
    return FreqDict2

def QualityNormal(fastqFile):
    from Bio import SeqIO
    import numpy as np
    import math
    ScoreList = np.array([])
    for record in SeqIO.parse(fastqFile, "fastq"):
        ScoreList = np.append(ScoreList,record.letter_annotations["phred_quality"])
    Error = 10**-(ScoreList/10)
    std = np.std(Error)
    mean = np.mean(Error)
    
    return mean

def positions(topN, lowN, freqDict):
    posList = []
    for key,value in freqDict.items():
        if value["N"]/sum(value.values()) <= lowN or value["N"]/sum(value.values()) >= topN:
            print(value["N"]/sum(value.values()))
            posList.append(key)
            print(freqDict[key])
    
    return posList


def VIPmain():
    from Clustering import plotit, fromTSV, Frequency, main
    import numpy as np
    
    tsv = sys.argv[1]
    fastqfile = sys.argv[2]
    high = QualityNormal(fastqfile)
    
    if len(sys.argv) > 3:
        assembly = sys.argv[3]
    else:
        assembly = False
    
    if assembly:
        FreqDict = main(assembly,fastqfile)
    else:
        AlignDict = fromTSV(tsv)
        FreqDict = Frequency(AlignDict)
    
    VIPList = VIPPos(FreqDict)
    FreqDict2 = NReplace(VIPList, fastqfile,tsv)
    NMean = plotit(FreqDict2)
    topN = NMean * (1+high)**3 
    lowN = NMean * (1-high)**3
    
    positionsList = (positions(topN,lowN, FreqDict2))
    
    
    return

VIPmain()