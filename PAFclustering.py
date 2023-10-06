#!/usr/bin/env python3
import sys
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram, cut_tree
import matplotlib.pyplot as plt




def PAFtoDM(paffilepath):
    # Initialize an empty dictionary to store read distances
    read_distances = {}

    with open(paffilepath, 'r') as f:
        for line in f:
                fields = line.strip().split('\t')
                read1, read2 = fields[0], fields[5]
                distance = int(fields[10]) / int(fields[9])
                
                # Add read distances to the dictionary
                if read1 not in read_distances:
                    read_distances[read1] = {}
                if read2 not in read_distances:
                    read_distances[read2] = {}
                
                if read1 == read2:
                    read_distances[read1][read2] = 0
                    read_distances[read2][read1] = 0
                else:
                    read_distances[read1][read2] = distance
                    read_distances[read2][read1] = distance
                
    read_distances = dict(sorted(read_distances.items()))
    # Create a DataFrame from the dictionary
    distance_matrix = pd.DataFrame.from_dict(read_distances, orient='index')
    print(distance_matrix)
   
    distance_matrix.fillna(float(10), inplace=True)
    
    print(distance_matrix)
    
    
    return distance_matrix



distance_matrix = PAFtoDM(sys.argv[1])
print(distance_matrix.shape)
linkage_matrix = linkage(distance_matrix, method='complete', metric='euclidean')
print(linkage_matrix)
cut = cut_tree(linkage_matrix, n_clusters = 3)
print(cut)
#print(distance_matrix)









