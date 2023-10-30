#!/usr/bin/env python3
import sys
import argparse
from collections import Counter
from Bio import SeqIO

def readBlast(blastfile):
    # Reads a blast results file and creates a dictionary with the entries 
    # with the highest alignmentscore for each contig that was blasted
    print("Reading blast results and saving to dict")
    blastdict = dict()
    
    with open(blastfile,'r') as file:
        
        seqcount = 0

        for line in file:
            seqcount += 1
            

            if not line.startswith("#"):
                line = line.split()
                readid = (line[0])

                if line[1].startswith("HLA"):
                    allele = line[1].split(":")[1].strip()
                    if readid not in blastdict:
                        blastdict[readid] = []
                    blastdict[readid].append(allele)
    print(seqcount)
            
    return blastdict

def readAllelelist(allelefile):
    alleledict = dict()
    with open(allelefile,'r') as file:
        for line in file:
            if not line.startswith("#"):
                line = line.split(",")
                alleleid = line[0]
                allele = line[1]
                alleledict[alleleid] = allele.strip()
    print(len(alleledict))
    return alleledict

def readHLAgen(hlanomfile):
    print("Reading HLA gene nomenclature G group file and returns it as a dict")
    # Reads the HLA gene nomenclature G group file and returns it as a dict
    hlanomdict = dict()
    HLAg = ""
    with open(hlanomfile,'r') as file:
        for line in file:
            if not line.startswith("#"):
                linesplit = line.strip().strip(";").split(";")
                HLAallele = linesplit[0].strip("*")
                if len(linesplit) > 2: #Checks if the G-group contains more than one entry
                    HLAlist = []
                    for gene in linesplit[1].split("/"):
                        HLAlist.append(gene)
                    HLAg = linesplit[2]
                    for gene in HLAlist:
                        hlanomdict[HLAallele,gene] = str(HLAallele + "*" + HLAg)
                else:
                    HLA = linesplit[1]
                    # If the line in the file only contains one entry, the HLA key is the same as its value since there is no G group for those genes.
                    hlanomdict[HLAallele,HLA] = str(HLAallele + "*" + HLA)
    
    return hlanomdict


def idtoallele(blastdict,alleledict,hlanomdict):
    ggroupdict = dict()
    finalggroup = dict()
    for readid in blastdict:
        for gene in blastdict[readid]:
            allele = tuple(alleledict[gene].split("*"))
            
            if allele[0] in ["A","B","C","DRB1","DQA1","DQB1","DPB1"]:
                
                if readid not in ggroupdict:
                    ggroupdict[readid] = []
                
                #print(hlanomdict)
                hlag = hlanomdict[allele]
                #print(hlag)
                #print(readid)
                ggroupdict[readid].append(hlag)
    #print(ggroupdict)
    for readid in ggroupdict:
        
        occurrence = Counter(ggroupdict[readid])
        
        
        finalggroup[readid] = occurrence
    #print(finalggroup)
    return ggroupdict

def genotyper(finalggroup):
    #print(finalggroup)
    allelecount = dict()
    for key in finalggroup:
        if finalggroup[key][0] not in allelecount:

            allelecount[finalggroup[key][0]] = 1
        else:
            allelecount[finalggroup[key][0]] += 1
    #print(allelecount)
    #print(sorted(allelecount.items(), key=lambda x:x[1]))
    count = dict()
    for i in allelecount.items():
        print(i)
        al = i[0].split(":")[0] #+ ":" + i[0].split(":")[1]
        if al not in count:
            
            count[al] = 1 * i[1]
        else:
            count[al] += 1 * i[1]
    cluster1, cluster2 = sorted(count.items(), key=lambda x:x[1])[-1], sorted(count.items(), key=lambda x:x[1])[-2]
    print(cluster1, cluster2)

    cluster1reads = set()
    cluster2reads = set()
    for key in finalggroup: 
        if finalggroup[key][0].split(":")[0] + ":" + finalggroup[key][0].split(":")[1] in cluster1:
            cluster1reads.add(key)
        if finalggroup[key][0].split(":")[0]  + ":" + finalggroup[key][0].split(":")[1] in cluster2:
            cluster2reads.add(key)
    print(len(cluster1reads))
    print(len(cluster2reads))
    print(cluster2reads)
    return cluster1reads, cluster2reads

    
def writeReads(readsets,fastqfile,outfile1,outfile2):

    if len(readsets) > 1:
        file1 = open(outfile1,'w')
        file2 = open(outfile2,'w')

        for read in SeqIO.parse(fastqfile, "fastq"):
            if read.id in readsets[0]:
                file1.write(read.format('fastq'))
            elif read.id in readsets[1]:
                file2.write(read.format('fastq'))
        
        file1.close()
        file2.close()   
        return
    else:
        file = open(outfile1,'w')
        for read in SeqIO.parse(fastqfile, "fastq"):
            file.write(read.format('fastq'))
        file.close()


class CustomParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def arguments():
    parser = CustomParser(description='Blastcluster.')
    parser.add_argument('--blastfile', type = str, required= True)
    parser.add_argument('--allelelist', type = str, required= True)
    parser.add_argument('--hlagnom', type = str, required= True)
    parser.add_argument('--fastq', type = str, required= True)
    parser.add_argument('--out1', type = str, required= True)
    parser.add_argument('--out2', type = str, required= True)

    return parser.parse_args()

def main():
    args = arguments()
    blastdict = readBlast(args.blastfile)
    alleledict = readAllelelist(args.allelelist)
    ggroupdict = readHLAgen(args.hlagnom)
    #matches = HLAmatcher(blastdict,ggroupdict)
    ggroups = (idtoallele(blastdict,alleledict,ggroupdict))
    clusters = genotyper(ggroups)
    writeReads(clusters, args.fastq, args.out1, args.out2)
main()