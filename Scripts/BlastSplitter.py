#!/usr/bin/env python3
import sys, os
import argparse
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
                readid = line[0]

                if line[1].startswith("HLA"):
                    allele = line[1].split(":")[1].strip()
                    if readid not in blastdict:
                        blastdict[readid] = []
                    blastdict[readid].append(allele)
                elif line[1] != "-":
                    allele = line[1].strip()
                    if readid not in blastdict:
                        blastdict[readid] = []
                    blastdict[readid].append(allele)
                
    #print(seqcount)
            
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
    #print(len(alleledict))
    return alleledict


def idtoallele(blastdict,alleledict):
    
    ggroupdict = dict()
    for readid in blastdict:
        for gene in blastdict[readid]:
            if not gene.startswith("HLA"):
                gene = f'HLA{gene}'
            try:
                allele = tuple(alleledict[gene].split("*"))
                
                if allele[0] in ["A","B","C","DRB1","DQA1","DQB1","DPB1"]:
                    
                    if readid not in ggroupdict:
                        ggroupdict[readid] = set()
                
                    ggroupdict[readid].add(allele[0])
            except KeyError as error:
                pass

    return ggroupdict

    
def writeReads(ggroupdict,fastqfile,ofA,ofB,ofC,ofDRB1,ofDQA1,ofDQB1,ofDPB1):
    print("Writing reads")
    if len(ggroupdict) > 0:
        ofA,ofBo,ofCo,ofDRB1o,ofDQA1o,ofDQB1o,ofDPB1o = open(ofA,'a'), open(ofB,'a'), open(ofC,'a'), open(ofDRB1,'a'), open(ofDQA1,'a'), open(ofDQB1,'a'), open(ofDPB1,'a')
        errorcount = 0
        successcount = 0
        #print(ggroupdict)
        for read in SeqIO.parse(fastqfile, "fastq"):
                if read.id in ggroupdict:
                    
                    ggroup = ggroupdict[read.id]
                    if ggroup == "A":
                        ofA.write(read.format('fastq'))
                    elif ggroup == "B":
                        ofB.write(read.format('fastq')) 
                    elif ggroup == "C":
                        ofC.write(read.format('fastq'))
                    elif ggroup == "DRB1":
                        ofDRB1.write(read.format('fastq'))
                    elif ggroup == "DQA1":
                        ofDQA1.write(read.format('fastq'))
                    elif ggroup == "DQB1":
                        ofDQB1.write(read.format('fastq'))
                    elif ggroup == "DPB1":
                        ofDPB1.write(read.format('fastq'))
                    else:
                        pass
                    successcount += 1
                else:
                    errorcount += 1
                    
                    
        print("Error count", errorcount)
        print("Success count", successcount)
        
        ofA.close(),ofBo.close(),ofCo.close(),ofDRB1o.close(),ofDQA1o.close(),ofDQB1o.close(),ofDPB1o.close()
    else: 
        print("No reads matched ggroup")
        return
    


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
    parser.add_argument('--oA', type = str, required= True)
    parser.add_argument('--oB', type = str, required= True)
    parser.add_argument('--oC', type = str, required= True)
    parser.add_argument('--oDRB1', type = str, required= True)
    parser.add_argument('--oDQA1', type = str, required= True)
    parser.add_argument('--oDQB1', type = str, required= True)
    parser.add_argument('--oDPB1', type = str, required= True)

    return parser.parse_args()

def main():
    args = arguments()
    print("Starting Blast Primer Splitting")
    blastdict = readBlast(args.blastfile)
    alleledict = readAllelelist(args.allelelist)
    #ggroupdict = readHLAgen(args.hlagnom)
    #matches = HLAmatcher(blastdict,ggroupdict)
    print("Blast Primer Splitting...")
    ggroups = (idtoallele(blastdict,alleledict))
   
    writeReads(ggroups, args.fastq, args.oA, args.oB,args.oC,args.oDRB1, args.oDQA1, args.oDQB1, args.oDPB1)
    print("Finished writing reads")
main()