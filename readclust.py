#!/usr/bin/env python3
import sys
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
#from scipy.cluster.hierarchy import linkage, dendrogram, cut_tree
#from scipy.spatial.distance import squareform


from collections import Counter
#from scipy.spatial.distance import squareform

class CustomParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def align(read1,read2):
    from Bio import pairwise2
    forwardAlign = pairwise2.align.globalms(read1, read2, 10, -5, -20, -20, score_only = True)
    revAlign = pairwise2.align.globalms(read1, read2.reverse_complement(), 10, -5, -20, -20, score_only = True)
    #forwardAlign = pairwisealigner.score(read1,read2,strand = '+')
    #revAlign = pairwisealigner.score(read1, read2.reverse_complement(), strand = '+')
    length = min(len(read1),len(read2))
    #rint(length)
    if forwardAlign >= revAlign:
        return 1/(forwardAlign/length)
       
    else:
        return 1/(revAlign/length)
    
#Create a python script that iterates over a fastq file and splits reads according to some criterion based on their name

def pairwiseAlignment(fastqfile,cutoff):
    
    
    read_distances = {}

    for read1 in SeqIO.parse(fastqfile, "fastq"):
        #print("Aligning seq:", read1)
        
        if len(read1) > cutoff:
        
            for read2 in SeqIO.parse(fastqfile, "fastq"):
                # Add read distances to the dictionary
                if read1.id not in read_distances:
                    read_distances[read1.id] = {}
                if read2.id not in read_distances:
                    read_distances[read2.id] = {}
                
                if len(read2) > cutoff:
                    
                    if read1.id == read2.id:
                        read_distances[read1.id][read2.id] = 0
                    else:
                        alignmentscore = align(read1,read2)
                        read_distances[read1.id][read2.id] = alignmentscore
                        read_distances[read2.id][read1.id] = alignmentscore
                        
   # Create a DataFrame from the dictionary
    distance_matrix = pd.DataFrame.from_dict(read_distances, orient='index')
    distance_matrix = distance_matrix.sort_index(axis=0).sort_index(axis=1)
    
    return distance_matrix

# def hierachichal(distancematrix):
#     linkage_matrix = linkage(distancematrix, method='complete', metric='euclidean')
#     #print(linkage_matrix)
#     cut = cut_tree(linkage_matrix, n_clusters = 2)
#     #print(cut)
#     return cut

def Agglomerative(distancematrix):
    from sklearn.cluster import AgglomerativeClustering
    distance_np = (distancematrix.to_numpy())
    clustering = AgglomerativeClustering().fit_predict(distance_np)
    print(sum(clustering), len(clustering)-sum(clustering))
    return clustering

def hdbscan_cluster(distance_matrix, min_cluster_size, min_samples = None):
    from sklearn.cluster import HDBSCAN
    import matplotlib.pyplot as plt
    import networkx as nx
    distance_matrix = (distance_matrix.to_numpy())
    
    G = nx.from_numpy_array(distance_matrix)
    print(G.edges(data=True))
    nx.draw(G, with_labels=True)
    plt.savefig('plotgraph.png', dpi=300, bbox_inches='tight')
    #print(G[1][1])

    norm = np.linalg.norm(distance_matrix, 1)
    distance_matrix = distance_matrix/norm


    print(distance_matrix)
    
    cluster = HDBSCAN(
        min_cluster_size = int(distance_matrix.shape[0] / 3),
        min_samples = int(distance_matrix.shape[0] / 3),
        metric='precomputed',
        cluster_selection_epsilon = 0.000000000000000,
        allow_single_cluster = True
    )
    
    # Perform HDBSCAN clustering
    clusterer = cluster.fit(distance_matrix)
    print(clusterer.labels_)   
    plt.imshow(distance_matrix, cmap='hot', interpolation='nearest')
    plt.savefig('distance_matrix.png')
    return clusterer.labels_



# From blast cluster
def readBlast(blastfile):
    # Reads a blast results file and creates a dictionary with the entries 
    # with the highest alignmentscore for each contig that was blasted
   
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


def labelingDistanceMatrix(distance_matrix:pd,alleledict,blastdict):
    print(distance_matrix)
    print(blastdict)
    for readname, _ in distance_matrix.iterrows():
        #print(col)
        #print(blastdict[readname][0])
        #print(readname.it))
        distance_matrix.at[readname,"Label"] = alleledict[blastdict[readname][0]]
    return distance_matrix

def divideReads(clustering,distancematrix):
    cluster0 = set()
    cluster1 = set()

    # Use Counter to count occurrences
    counts = Counter(clustering)
    # Use most_common to get the two most common numbers
    
    most_common = counts.most_common(2)
    if (len(most_common) > 1):
        print(most_common[1][1]/(most_common[0][1]+most_common[1][1]))
        print(clustering)
        for i in range(len(clustering)):
            if clustering[i] == most_common[0][0]:
                cluster0.add(distancematrix.index[i])
            elif clustering[i] == most_common[1][0]:
                cluster1.add(distancematrix.index[i])
        print("Sizes:", len(cluster0), len(cluster1))
        return cluster0, cluster1
    else:
        print("Found only one cluster.")
        return set(clustering)
        

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

def arguments():
    parser = CustomParser(description='Pairwise alignment (AVA).')
    parser.add_argument('--fastq', type = str, required= True)
    parser.add_argument('--csv', type = str, required = True)
    parser.add_argument('--out1', type = str)
    parser.add_argument('--out2', type = str)
    parser.add_argument('--blastresult', type = str)
    parser.add_argument('--allelelist', type = str)
    parser.add_argument('--mode', type = str, choices=['cm', 'fcsv', 'almat'], required = True)

    
    return parser.parse_args()


def continousmode():
    args = arguments()
    distancematrix = (pairwiseAlignment(args.fastq, cutoff  = 3000))
    distancematrix.to_csv(args.csv)
    aggclustering = Agglomerative(distancematrix)
    print(aggclustering)
    hdbscan_cluster(distancematrix, min_cluster_size = 10, min_samples= None)
    readsets = divideReads(aggclustering,distancematrix)

    
    writeReads(readsets,args.fastq,args.out1,args.out2)
    
def aligntomatrix():
    args = arguments()
    distancematrix = (pairwiseAlignment(args.fastq, cutoff  = 3000))
    distancematrix.to_csv(args.csv)

def runfromCSV():
    args = arguments()
    distancematrix = pd.read_csv(args.csv,header=0, index_col=0)
    #print(distancematrix)
    #aggclustering = Agglomerative(distancematrix)
    blastdict = readBlast(args.blastresult)
    alleledict = readAllelelist(args.allelelist)
    distancematrix_labelled = labelingDistanceMatrix(distancematrix,alleledict,blastdict)

    hdbclustering = hdbscan_cluster(distancematrix_labelled, min_cluster_size = 40, min_samples= None)
   # print(aggclustering)
    

    readsets = divideReads(hdbclustering,distancematrix)
    
    writeReads(readsets,args.fastq,args.out1,args.out2)
    #print(distancematrix)
    

def main():
    args = arguments()
    if args.mode == "cm":
        #continousmode()
        aligntomatrix()
        runfromCSV()
    elif args.mode == "fcsv":
        runfromCSV()
    elif args.mode == "almat":
        aligntomatrix()
main()