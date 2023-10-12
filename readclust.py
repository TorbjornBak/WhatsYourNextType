#!/usr/bin/env python3
import sys
import argparse
import pandas as pd
from Bio import SeqIO
from Bio import pairwise2
import numpy as np

from scipy.cluster.hierarchy import linkage, dendrogram, cut_tree
from sklearn.cluster import AgglomerativeClustering

from progress.bar import Bar

class CustomParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def align(read1,read2):
    forwardAlign = pairwise2.align.globalxx(read1, read2, score_only = True)
    revAlign = pairwise2.align.globalxx(read1, read2.reverse_complement(), score_only = True)
    #forwardAlign = pairwisealigner.score(read1,read2,strand = '+')
    #revAlign = pairwisealigner.score(read1, read2.reverse_complement(), strand = '+')
    
    if forwardAlign >= revAlign:
        return 1/forwardAlign
       
    else:
        return 1/revAlign
    
#Create a python script that iterates over a fastq file and splits reads according to some criterion based on their name

def pairwiseAlignment(fastqfile,cutoff):
    read_distances = {}
    fastq_dict = SeqIO.index(fastqfile, "fastq")
    bar = Bar('Aligning sequences', max=len(fastq_dict))
   
    for read1 in SeqIO.parse(fastqfile, "fastq"):
        #print("Aligning seq:", read1)
        bar.next()
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
                        read_distances[read2.id][read1.id] = 0
                    else:
                        alignmentscore = align(read1,read2)
                        read_distances[read1.id][read2.id] = alignmentscore
                        read_distances[read2.id][read1.id] = alignmentscore
    
    bar.finish()
    read_distances = dict(sorted(read_distances.items()))
    # Create a DataFrame from the dictionary
    distance_matrix = pd.DataFrame.from_dict(read_distances, orient='index')
    
    return distance_matrix

def hierachichal(distancematrix):
    linkage_matrix = linkage(distancematrix, method='complete', metric='euclidean')
    #print(linkage_matrix)
    cut = cut_tree(linkage_matrix, n_clusters = 2)
    #print(cut)
    return cut

def Agglomerative(distancematrix):
    clustering = AgglomerativeClustering().fit_predict(distancematrix)
    print(sum(clustering), len(clustering)-sum(clustering))
    return clustering

def hdbscan_cluster(distance_matrix, min_cluster_size, min_samples = None):
    cluster = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        metric='precomputed'  # Use precomputed distances
    
    )
    
    # Perform HDBSCAN clustering
    clusterer = cluster.fit(distance_matrix)

    distance_matrix['ClusterLabel'] = clusterer.labels_

    # Calculate the sizes of each cluster
    cluster_sizes = distance_matrix['ClusterLabel'].value_counts()

    # Sort the clusters by size in descending order
    sorted_clusters = cluster_sizes.sort_values(ascending=False)

    # Select the two largest clusters
    if len(sorted_clusters) >= 2:
        largest_cluster = sorted_clusters.index[0]
        second_largest_cluster = sorted_clusters.index[1]
        print("Largest Cluster Label:", largest_cluster)
        print("Size of Largest Cluster:", sorted_clusters.iloc[0])
        print("Second Largest Cluster Label:", second_largest_cluster)
        print("Size of Second Largest Cluster:", sorted_clusters.iloc[1])
    else:
        print("There are less than two clusters.")

    # If you want to get the data points in the largest clusters, you can use:
    largest_cluster_points = distance_matrix[distance_matrix['ClusterLabel'] == largest_cluster]
    second_largest_cluster_points = distance_matrix[distance_matrix['ClusterLabel'] == second_largest_cluster]
    print(largest_cluster_points)
    print(second_largest_cluster_points)
        

def divideReads(clustering,distancematrix):
    cluster0 = set()
    cluster1 = set()
    for i in range(len(clustering)):
        if clustering[i] == 0:
            cluster0.add(distancematrix.index[i])
        elif clustering[i] == 1:
            cluster1.add(distancematrix.index[i])
    return cluster0, cluster1
        
def writeReads(readsets,fastqfile,outfile1,outfile2):

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

def arguments():
    parser = CustomParser(description='Pairwise alignment (AVA).')
    parser.add_argument('--fastq', type = str, required= True)
    parser.add_argument('--csv', type = str, required = True)
    parser.add_argument('--out1', type = str)
    parser.add_argument('--out2', type = str)
    parser.add_argument('--mode', type = str, choices=['cm', 'fcsv', 'almat'], required = True)
    
    
   
    return parser.parse_args()

def continousmode():
    args = arguments()
    distancematrix = (pairwiseAlignment(args.fastq, cutoff  = 2000))
    distancematrix.to_csv(args.csv)
    aggclustering = Agglomerative(distancematrix)
    print(aggclustering)

    readsets = divideReads(aggclustering,distancematrix)
    print("Sizes:", len(readsets[0]), len(readsets[1]))
    writeReads(readsets,args.fastq,args.out1,args.out2)
    
def aligntomatrix():
    args = arguments()
    distancematrix = (pairwiseAlignment(args.fastq, cutoff  = 2000))
    distancematrix.to_csv(args.csv)

def runfromCSV():
    args = arguments()
    distancematrix = pd.read_csv(args.csv,header=0, index_col=0)
    aggclustering = Agglomerative(distancematrix)
    print(aggclustering)

    readsets = divideReads(aggclustering,distancematrix)
    print("Sizes:", len(readsets[0]), len(readsets[1]))
    writeReads(readsets,args.fastq,args.out1,args.out2)

def main():
    args = arguments()
    if args.mode == "cm":
        continousmode()
    elif args.mode == "fcsv":
        runfromCSV()
    elif args.mode == "almat":
        aligntomatrix()
main()