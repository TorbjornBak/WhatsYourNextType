#!/usr/bin/env python3
import sys
import argparse
import subprocess
import os.path
import pandas as pd
import glob

def readBlast(blastfile):
    # Reads a blast results file and creates a dictionary with the entries 
    # with the highest alignmentscore for each contig that was blasted
    blastdict = dict()
    Query = ""
    HLAvariant = ""
    AlignmentLength = 0
    AlignmentScore = 0
    with open(blastfile,'r') as file:
        HLAflag = False
        seqcount = 0
  
        for line in file:
            
            if line.startswith("Query="):
                seqcount += 1
                
                Query = str(seqcount)
                HLAflag = False
                AlignmentLength = None
                #lowestLength = None
                AlignmentScore = None
                prevAlignmentScore = None
                       
            elif line.startswith("HLA") and HLAflag == False and (AlignmentScore is None or int(line.split()[4])*1.001 >= int(prevAlignmentScore)):
                linesplit = line.split()
                HLAvariant = [linesplit[1].split("*")[0], linesplit[1].split("*")[1]]
                    
                    
                AlignmentScore = int(linesplit[4])
                alignmentScoreCutoff = 3000
                
                if AlignmentScore > alignmentScoreCutoff: 
                    AlignmentLength = int(linesplit[2])
                    
                    if Query not in blastdict:
                        blastdict[Query] = {}

                    if (prevAlignmentScore == None or int(AlignmentScore)*1.001 >= int(prevAlignmentScore)) and len(blastdict[Query]) < 3:
                        HLAflag = False
                        
                    else: 
                        HLAflag = True
                        
                    variantNr = f'V{len(blastdict[Query])}'
                    blastdict[Query][variantNr] = {'HLAvariant':HLAvariant,'AlignmentLength':AlignmentLength, 'AlignmentScore':AlignmentScore}  
                    
                # if lowestLength is None or AlignmentLength < lowestLength:
                    #    likelyVariant = variantNr
                    #   lowestLength = AlignmentLength
                    #  blastdict[Query]["Variant"] = blastdict[Query][likelyVariant]
                    
                
                prevAlignmentScore = AlignmentScore  

            
    
   
   
    return blastdict

def readHAPDUP(hapdupLog):
    #Reads HAPDUP / margin log and returns a dict with the information regarding haplotype distribution
    import re
    HapDupDict = dict()
    with open(hapdupLog, "r") as file:
        for line in file:
            match = re.search(r"(.*)_bin.*", line)
            if match.group(1)[0:3] == "HLA":
                gene = match.group(1)[3]
            else:
                gene = match.group(1)
            HapDupDict[gene] = line.strip().split(", ")[1:]
    return HapDupDict




def readHLAgen(hlanomfile):
    #print("Reading HLA gene nomenclature G group file and returns it as a dict")
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
                        hlanomdict[HLAallele,gene] = [HLAallele, HLAg]
                else:
                    HLA = linesplit[1]
                    # If the line in the file only contains one entry, the HLA key is the same as its value since there is no G group for those genes.
                    hlanomdict[HLAallele,HLA] = [HLAallele, HLA]
                    
    return hlanomdict

def HLAmatcher(blastdict,hlanomdict, allelelist):
    #print("Matching HLA genes to the HLA gene nomenclature G group dictionary")
    #Matches the blast results to the respective g groups from the hla nom dict and saves the genes to a list
    #genes = list()
    #print(blastdict)
    df = pd.DataFrame()
    grouplist = list()
    typelist = list() 
    exactTypeList = list()
    alternativeTypesList = list()
    scoreList = list()
    for gene in blastdict:
        if blastdict[gene]['V0']['HLAvariant'][0] in allelelist:
           
            #blastdict -> querynr -> variant -> {'HLAvariant':HLAvariant,'AlignmentLength':AlignmentLength, 'AlignmentScore':AlignmentScore}
            hlatype = blastdict[gene]['V0']['HLAvariant']
            translatedHlaType = hlanomdict[hlatype[0],hlatype[1]]
            score = blastdict[gene]['V0']['AlignmentScore']
            altList = list()
            for HLAtype in blastdict[gene]:
                
                altType = blastdict[gene][HLAtype]["HLAvariant"]
                
                
                if HLAtype != 'V0' and hlanomdict[altType[0],altType[1]] != translatedHlaType:
                    #print(altType)
                    altList.append(hlanomdict[altType[0],altType[1]])
           
            alternativeTypesList.append(altList)
            grouplist.append(translatedHlaType[0])
            typelist.append(translatedHlaType[1])
            exactTypeList.append("*".join(hlatype))
            scoreList.append(score)
            
            
    
    df["HLAgroup"] = grouplist
    df["HLAnr"] = typelist
    df["HLAtypeG"] = [f'{g}*{t}' for g, t in zip(grouplist, typelist)]
    df["ExactHLAtype"] = exactTypeList
    df["AlternativeTypes"] = alternativeTypesList
    df["Score"] = scoreList
    
    lastGroup = None
    lastList = list()
    HNrList = list()
    for group in df["HLAgroup"]:
        if lastGroup != group:
            
            lastGroup = group
            lastList = [group]
            HNrList.append(f"H{len(lastList)}")
        else:
            lastList.append(group)
            HNrList.append(f"H{len(lastList)}")
            lastGroup = group
            
    df["HNr"] = HNrList
    #print(df)

    return df

#def missingAlleles(genes,allelelist):
    print("Checking if any alleles are missing.")
    # Checking if there is at least one of each HLA gene 
    # and checks that there is not more than two of each
    genecount = [0] * len(allelelist) # Initalize empty counter
    for gene in genes:
        index = allelelist.index(gene[0][0])
        genecount[index] += 1
    errorcheck = False
    for i in range(len(genecount)):
        if genecount[i] > 2:
            #print("Error: More than two alleles of gene", allelelist[i], "was found. Do not use these results!")
            errorcheck = True
        elif genecount[i] == 0:
            print("Error: No alleles of gene", allelelist[i], "was found. Do not use these results!")
            errorcheck = True
    if errorcheck == False:
        print("There was at least one of each gene and not more than two of each. Passed.")
    print(genecount)  
    return allelelist, genecount

def missingAllelesPD(HLAdf,allelelist):
    #print("Checking if any alleles are missing.")
    # Checking if there is at least one of each HLA gene 
    # and checks that there is not more than two of each
    genecount = [0] * len(allelelist) # Initalize empty counter
    #print(HLAdf["HLAgroup"])
    for gene in list(HLAdf["HLAgroup"]):
        index = allelelist.index(gene)
        genecount[index] += 1
    errorcheck = False
    for i in range(len(genecount)):
        if genecount[i] > 2:
            #print("Error: More than two alleles of gene", allelelist[i], "was found. Do not use these results!")
            errorcheck = True
        elif genecount[i] == 0:
            #print("Error: No alleles of gene", allelelist[i], "was found. Do not use these results!")
            errorcheck = True
    if errorcheck == False:
        print("Success: there were at least one of each gene and not more than two of each.")
    else:
        print("Error")
    print(genecount)  
    
    return allelelist, genecount

def benchmarking(correcttype,predictedtypes):
    #NOT IMPLEMENTED
    pass


# def printGenes(genes, marginLog):
#     printString = []
#     tempGene = ""
#     for gene in genes:
#         if gene not in printString:
#             printString.append(gene)
#     HapDupDict = readHAPDUP(marginLog)
#     print("G-group \t\tAllelebalance")
#     iterator = 0
#     for gene in sorted(printString):
#         # if gene != tempGene:
#         tempGene = gene[0]
#         if tempGene[0] in HapDupDict.keys():
#             print(tempGene,end = "\t")
#             print("H"+str(iterator+1), HapDupDict[tempGene[0]][iterator])
#             iterator += 1
#         else:
#             print(gene[0], "UNPHASED")
#         if tempGene[0] in HapDupDict.keys():
#             if iterator == len(HapDupDict[tempGene[0]]):
#                 iterator= 0

#     return

# def printAllGenes(genes):
#     print("Showing all found genes: ")
#     for gene in genes:
#         print(gene)
#     return

#def printToFile(genes, outputpath, marginLog):
    #Saving to file
    printString = []
    printDict = dict() 
    tempGene = []
    for gene in genes:
        if gene[0] != tempGene:
            tempGene = gene[0]
            if gene[0][0] not in printDict.keys():
                printDict[gene[0][0]] = dict()
                printDict[gene[0][0]]["G-group"] = list()
            printDict[gene[0][0]]["G-group"].append(gene[0])
            if gene not in printString:
                printString.append(gene)
    outputfile = open(outputpath,'w')
    HapDupDict = readHAPDUP(marginLog)
    
    for genes,gGroup in printDict.items():
        for groups in gGroup["G-group"]:
            print(groups, file = outputfile, end = "\t")
        if genes in HapDupDict.keys():
            for i in range(len(HapDupDict[genes])):
                NumList = [1,2,0]
                print("H"+str(NumList[i]) +":"+HapDupDict[genes][i], file = outputfile, end = "\t")
        print("", file = outputfile)


    # for gene in sorted(printString):
    #     if gene[0][0] != tempGene:
    #         tempGene = gene[0][0]"
    #         print(gene,"\t",file = outputfile, end ="\t")
    #         if tempGene in HapDupDict.keys():
    #             for i in range(len(HapDupDict[tempGene])):
    #                 NumList = [1,2,0]
    #                 print("H"+str(NumList[i]) +":"+HapDupDict[tempGene][i], file = outputfile, end = "\t")
    #         print("", file = outputfile)
        # for gene in sorted(printString):
        #     
    outputfile.close()
    return

def addHaplotypeInfo(haplotypediv, df):
    haplotypedivList = list()
    unmappedreadsList = list()
    for group, hnr in zip(df["HLAgroup"], df["HNr"]):
       
        if group in haplotypediv and len(haplotypediv[group]) > 1:
            if hnr == "H1":
                haplotypedivList.append(haplotypediv[group][0])
            else:
                haplotypedivList.append(haplotypediv[group][1])
            unmappedreadsList.append(haplotypediv[group][2])
        else:
            haplotypedivList.append("Unphased")
            unmappedreadsList.append("-")
    
    df["HapDupPhaseCoverage"] = haplotypedivList
    df["UnmappedReadCoverage"] = unmappedreadsList

    return df

def sortDataFramebyList(df, sorted_list):
    df['sort_index'] = df['HLAgroup'].map({val: i for i, val in enumerate(sorted_list)})
    df_sorted = df.sort_values(by='sort_index').drop('sort_index', axis=1)
    df_sorted = df_sorted.reset_index(drop=True)
    return df_sorted

class CustomParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def arguments():
    parser = CustomParser(description='Genotyper for HLA.')
    parser.add_argument('--blastfile', type = str, help = 'specify the blast results file')
    parser.add_argument('--hlagen', type = str, help = 'specify the HLA gen g group file')
    parser.add_argument('--reads', type = str)
    parser.add_argument('--hladatabase', type = str)
    parser.add_argument('--primerlist', type = str)
    parser.add_argument('--threads', type = int, default = 4)
    parser.add_argument('--output', type = str)
    parser.add_argument('--marginLog', default= None, type = str)
    parser.add_argument('--recursive', default = False, type = bool)
    
    
    return parser.parse_args()


def recursive_main(blastpaths, hlagen, marginLog):
    print(blastpaths)
    blastpaths = glob.iglob(blastpaths)
    print(blastpaths)
    for blastfile in blastpaths:
        blastdict = readBlast(blastfile)
        hlanomdict = readHLAgen(hlagen)
        allelelist = ["A","B","C","DRB1","DQA1","DQB1","DPB1"]
        resultDF = HLAmatcher(blastdict,hlanomdict,allelelist)

        missingAllelesPD(resultDF,allelelist)
        if marginLog is not None:
            haplotypediv = readHAPDUP(marginLog)
            resultDF = addHaplotypeInfo(haplotypediv, resultDF)
        
        resultDF = sortDataFramebyList(resultDF,allelelist)
        outputpath = os.path.splitext(blastfile)[0] + '_HLA_type.tsv'
        resultDF.to_csv(outputpath, sep = '\t')

def main():
    args = arguments()

    if args.recursive is True:
        return recursive_main(args.blastfile, args.hlagen, args.marginLog)
    
    blastdict = readBlast(args.blastfile)
    hlanomdict = readHLAgen(args.hlagen)
    allelelist = ["A","B","C","DRB1","DQA1","DQB1","DPB1"]
    resultDF = HLAmatcher(blastdict,hlanomdict,allelelist)

    missingAllelesPD(resultDF,allelelist)
    if args.marginLog is not None:
        haplotypediv = readHAPDUP(args.marginLog)
        resultDF = addHaplotypeInfo(haplotypediv, resultDF)
    
    resultDF = sortDataFramebyList(resultDF,allelelist)
    resultDF.to_csv(args.output, sep = '\t')
    print(resultDF)


main()

'''def testing():
    blastdict = readBlast(sys.argv[1])
    hlanomdict = readHLAgen(sys.argv[2])
    genes = HLAmatcher(blastdict,hlanomdict)
    allelelist = ["A","B","C","DRB1","DQA1","DQB1","DPB1"]
    missingAlleles(genes,allelelist)
    #print(blastdict)

    #print(hlanomdict)
    printGenes(genes)

    printAllGenes(genes)
testing()'''