#!/usr/bin/env python3
import sys
import argparse
import subprocess
import os.path


def assembleReads(readfiles,threads):
    print("Assembling reads")
    resultdir = "Assembly/"
    subprocess.run(["mkdir",resultdir, "-p"])
    for reads in readfiles:
        outputpath = resultdir + os.path.splitext(reads)[0]

        subprocess.run(['flye', "--nano-hq", reads, "--read-error", "0.05", "--threads", str(threads), "--out-dir", outputpath])
       # subprocess.run(['minimap2', '-ax', 'map-ont', '-t', str(threads), outputpath, reads, '|', 'samtools', 'sort', '-@', '4', '-m', '4G', '>', 'lr_mapping.bam'],shell=True, check=True)
      #  subprocess.run(['samtools', 'index', '-@', '4', 'lr_mapping.bam'])
     #   pwdir = f'{resultdir}:{resultdir}'
    #    assemblypath = outputpath +"/assembly.fasta"
   #     subprocess.run(['docker', 'run', '-v', pwdir, '-u', ""`id"", ""-u`:`id"", ""-g`"", "mkolmogo/hapdup:0.12", "\",  "hapdup", "--assembly", assemblypath, "--bam", "HD_DIR/lr_mapping.bam", "--out-dir", "HD_DIR/hapdup", "-t", threads, "--rtype", "hifi"])
        print("Finished assembly for:", reads)
   
    
    #for file in readfiles:
    #   subprocess.run(['rm',file])

    mergepattern = resultdir + "*/assembly.fasta"
    mergedfasta = "all_assemblies.fasta"
    mergecommand = f'cat {mergepattern} > {mergedfasta}'

    print("Merging files to one fasta:", mergedfasta)
    subprocess.run(mergecommand, shell=True, check=True)
    
    return mergedfasta
    

def blastSequences(fastafile,databasename,threads):
    print("Blasting sequences")
    outputpath = "blast" + os.path.splitext(fastafile)[0] + '.txt'
    subprocess.run(["blastn", "-query", fastafile, "-db", databasename, "-out", outputpath, "-num_alignments", "3", "-gapopen", "3", "-gapextend", "2", "-penalty", "-3", "-reward", "1", "-word_size", "11", "-num_threads", str(threads)])
    return outputpath


def readBlast(blastfile):
    # Reads a blast results file and creates a dictionary with the entries 
    # with the highest alignmentscore for each contig that was blasted
    print("Reading blast results and saving to dict")
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
                #Query = line.split()[1]
                Query = str(seqcount)
                HLAflag = False
                       
            elif line.startswith("HLA") and HLAflag == False:
                linesplit = line.split()
                HLAvariant = [linesplit[1].split("*")[0], linesplit[1].split("*")[1]]
                AlignmentLength = linesplit[2]
                AlignmentScore = linesplit[4]
                blastdict[Query] = [HLAvariant,AlignmentLength,AlignmentScore]
                HLAflag = True
                seqcount += 1
                
            
    return blastdict


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
                        hlanomdict[HLAallele,gene] = [HLAallele, HLAg]
                else:
                    HLA = linesplit[1]
                    # If the line in the file only contains one entry, the HLA key is the same as its value since there is no G group for those genes.
                    hlanomdict[HLAallele,HLA] = [HLAallele, HLA]
    return hlanomdict

def HLAmatcher(blastdict,hlanomdict):
    print("Matching HLA genes to the HLA gene nomenclature G group dictionary")
    #Matches the blast results to the respective g groups from the hla nom dict and saves the genes to a list
    
    genes = list()
    for gene in blastdict:
        genes.append(hlanomdict[blastdict[gene][0][0],blastdict[gene][0][1]])
    
    return sorted(genes)

def missingAlleles(genes,allelelist):
    print("Checking if any alleles are missing.")
    # Checking if there is at least one of each HLA gene 
    # and checks that there is not more than two of each
    genecount = [0] * len(allelelist) # Initalize empty counter
    for gene in genes:
        index = allelelist.index(gene[0])
        genecount[index] += 1
    errorcheck = False
    for i in range(len(genecount)):
        if genecount[i] > 2:
            print("Error: More than two alleles of gene", allelelist[i], "was found. Do not use these results!")
            errorcheck = True
        elif genecount[i] == 0:
            print("Error: No alleles of gene", allelelist[i], "was found. Do not use these results!")
            errorcheck = True
    if errorcheck == False:
        print("There was at least one of each gene and not more than two of each. Passed.")
    print(genecount)  
    return allelelist, genecount

def benchmarking(correcttype,predictedtypes):

    pass


def printGenes(genes):
    printString = []
    for gene in genes:
        if gene not in printString:
            printString.append(gene)
    for gene in sorted(printString):
        print(gene[0],"\t",gene[1])
    return

def printAllGenes(genes):
    print("Showing all found genes: ")
    for gene in genes:
        print(gene[0],"\t",gene[1])
    return

def printToFile(genes, outputpath):
    #Saving to file
    printString = []
    for gene in genes:
        if gene not in printString:
            printString.append(gene)
    outputfile = open(outputpath,'w')

    for gene in sorted(printString):
        print(gene[0],"\t",gene[1],file = outputfile)
    outputfile.close()
    return



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
    
    return parser.parse_args()

# def main():
#     args = arguments()
    
#     readfiles = PrimerSplitter.PrimerSplitter(args.primerlist,args.reads)
#     assembly = assembleReads(readfiles,args.threads)

#     blastfile = blastSequences(assembly,args.hladatabase,args.threads)
#     blastdict = readBlast(blastfile)
#     hlanomdict = readHLAgen(args.hlagen)
#     genes = HLAmatcher(blastdict,hlanomdict)
#     allelelist = ["A","B","C","DRB1","DQA1","DQB1","DPB1"]
#     missingAlleles(genes,allelelist)

#     printGenes(genes)

#     printAllGenes(genes)

# main()

def main():
    args = arguments()
    
    blastdict = readBlast(args.blastfile)
    hlanomdict = readHLAgen(args.hlagen)
    genes = HLAmatcher(blastdict,hlanomdict)
    allelelist = ["A","B","C","DRB1","DQA1","DQB1","DPB1"]
    
    missingAlleles(genes,allelelist)

    printGenes(genes)

    printAllGenes(genes)

    printToFile(genes,args.output)

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