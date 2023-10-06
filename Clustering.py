import numpy as np
import pandas as pd
import sys



file = sys.argv[1]
fileName = sys.argv[2]
fastQFileIn = sys.argv[3]

header = [ "query len","query start", "query end", "+/-", "target name", "target length", "target start", "target end","matching base", "gaps", "mapping quality","1","2","3","4","5"]
data = pd.read_csv(file, sep= "\t", names = header)


data["dist"] = data["matching base"]/ data["gaps"]

avgDist = np.mean(data["dist"])
Cluster1 = data.index[data["dist"] >= avgDist]
Cluster2 = data.index[data["dist"] < avgDist]


FastQFile = open(fastQFileIn,"r")
OutputFile1 = open(fileName+".cluster1.fastq","w")
OutputFile2 = open(fileName+".cluster2.fastq","w")

FastaFlag1 = False
FastaFlag2 = False
iterator = 0
for line in FastQFile:
    if line.split(" ")[0][1:] in Cluster1:
        FastaFlag1 = True
        iterator = 1
    if line.split(" ")[0][1:] in Cluster2:
        FastaFlag2 = True
        iterator = 1
    
    if FastaFlag1 == True:
        print(line, file=OutputFile1, end= "")
    if FastaFlag2 == True:
        print(line, file=OutputFile2, end="")
        
    if iterator == 4:
        FastaFlag1 = False
        FastaFlag2 = False
    iterator += 1
    