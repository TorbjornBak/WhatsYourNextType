
import multiprocessing


def reverseComplementaryprimer(primer):
    reverseStrand = ""
    for char in primer.upper():
        if char == "A":
            reverseStrand += "T"
        if char == "T":
            reverseStrand += "A"
        if char == "G":
            reverseStrand += "C"
        if char == "C":
            reverseStrand += "G"
            
    return reverseStrand[::-1]


def PrimerSplitter(PrimerList, FastqFile):


    
    import subprocess
    infile = open(PrimerList, "r")
    PrimerDict = dict()
    for line in infile:
        
        lineTMP = line.strip().split("\t")
        
        if "HLA" in lineTMP[0].strip().split("-"):

            if lineTMP[0].split("-")[0]+lineTMP[0].split("-")[1] in PrimerDict.keys():
                PrimerDict[lineTMP[0].split("-")[0]+lineTMP[0].split("-")[1]].append(lineTMP[1].upper())
            else:
                PrimerDict[lineTMP[0].split("-")[0]+lineTMP[0].split("-")[1]] = []
                PrimerDict[lineTMP[0].split("-")[0]+lineTMP[0].split("-")[1]].append(lineTMP[1].upper())

        if line[0] == "D":
            if lineTMP[0].split("-")[0] in PrimerDict.keys():
                
                PrimerDict[lineTMP[0].split("-")[0]].append(lineTMP[1].upper())
            else: 
                PrimerDict[lineTMP[0].split("-")[0]] = []
                PrimerDict[lineTMP[0].split("-")[0]].append(lineTMP[1].upper())
    infile.close()

    from Bio import SeqIO
    
    for keys1 in PrimerDict.keys():
        subprocess.run(["rm", keys1+"_bin.fastq", "-f"])
    recordadded = False
    printFlag = False
    FileList = set()
    for record in SeqIO.parse(FastqFile, "fastq"):
        for keys in PrimerDict.keys():
            for primer in PrimerDict[keys]:
                for i in range(4):
                    for b in range(4):
                        if record.seq[b:len(primer)+b-i] == (primer[b:-i] or primer[i:] or primer[0:-b]) or record.seq[0:-b] == primer[b:] or record.seq[0:len(primer)-b] == primer[b:]:
                            recordadded = True
                            writefile = open(workDir+"/Bins/"+keys+"_bin.fastq", "a")
                            FileList.add(keys+"_bin.fastq")
                            print("@"+record.id, file = writefile)
                            print(record.seq[15:], file = writefile)
                            print("+", file = writefile)
                            print("".join(chr(q + 33) for q in record.letter_annotations["phred_quality"])[15:], file = writefile)
                            break
                        
                        
                        revPrimer = reverseComplementaryprimer(primer)
                        #revPrimer = primer.reverse_complement()
                        if record.seq[len(record.seq)-len(primer)- b:len(record.seq)-b] == (revPrimer[0:] or revPrimer[0:-b] or revPrimer[b:]):
                            recordadded = True
                            writefile = open(workDir+"/Bins/"+keys+"_bin.fastq", "a")
                            FileList.add(keys+"_bin.fastq")
                            print("@"+record.id, file = writefile)
                            print(record.seq[:-15], file = writefile,)
                            print("+", file = writefile)
                            print("".join(chr(q + 33) for q in record.letter_annotations["phred_quality"])[:-15], file = writefile)
                            break
                            
                    if recordadded == True:
                        break
                if recordadded == True:
                        break
            if recordadded == True:
                    break
        if recordadded == False:
            writefile = open(workDir+"/Bins/excess_bin.fastq", "a")
            print(record.format("fastq"), file = writefile, end ="")       
        recordadded = False
    
    FileList = list(FileList)
    FileList = PrimerAligner(PrimerDict, FileList, workDir+"/Bins/excess_bin.fastq")
    #subprocess.run(["rm",workDir+"/Bins/excess_bin.fastq"])
    return(FileList)

def PrimerAligner(PrimerDict, FileList,fastq = "excess_bin.fastq",):
    alignmentcounter = 0
    from Bio import SeqIO
    from Bio import pairwise2
    Align1Flag = True
    for record in SeqIO.parse(fastq, "fastq"):
        maxAlignmentScore = 0
        maxKey = ""
        for keys in PrimerDict.keys():
            for primer in PrimerDict[keys]:
                revPrimer = reverseComplementaryprimer(primer)
                #revPrimer = primer.reverse_complement()
                revRecord = record.seq[len(record.seq)-len(primer)-4:]
                
                seq2 = record.seq[0:len(primer)+4]
                forwardAlign = pairwise2.align.globalxx(primer,seq2, score_only = True)
                revAlign = pairwise2.align.globalxx(revPrimer,revRecord, score_only = True)
                
                if forwardAlign > revAlign: 
                    if forwardAlign > maxAlignmentScore:
                        maxAlignmentScore = forwardAlign
                        maxKey = keys
                        Align1Flag = True

                else:
                    if revAlign > maxAlignmentScore:
                        maxAlignmentScore = revAlign
                        maxKey = keys
                        Align1Flag = False
                
                

        if maxAlignmentScore > 18:
            writefile = open(workDir+"/Bins/"+maxKey+"_bin.fastq","a")
            if Align1Flag == True:
                print("@"+record.id, file = writefile)
                print(record.seq[15:], file = writefile)
                print("+", file = writefile)
                print("".join(chr(q + 33) for q in record.letter_annotations["phred_quality"])[15:], file = writefile)
            else:
                print("@"+record.id, file = writefile)
                print(record.seq[:-15], file = writefile)
                print("+", file = writefile)    
                print("".join(chr(q + 33) for q in record.letter_annotations["phred_quality"])[:-15], file = writefile)
            
            writefile.close()
        else:
            writefile = open(workDir+"/Bins/excess2_bin.fastq","a")
            if Align1Flag == True:
                print("@"+record.id, file = writefile)
                print(record.seq[15:], file = writefile)
                print("+", file = writefile)
                print("".join(chr(q + 33) for q in record.letter_annotations["phred_quality"])[15:], file = writefile)
                
                
            else:
                print("@"+record.id, file = writefile)
                print(record.seq[:-15], file = writefile)
                print("+", file = writefile)    
                print("".join(chr(q + 33) for q in record.letter_annotations["phred_quality"])[:-15], file = writefile)

                
            writefile.close()
    return(FileList)
import sys

PrimerList = sys.argv[1]
FastQFile = sys.argv[2]
workDir = sys.argv[3]

PrimerSplitter(PrimerList, FastQFile)


        
        
    
