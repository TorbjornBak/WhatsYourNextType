
import pandas as pd
import sys
import argparse
#from readcllust import hierachichal
from scipy.cluster.hierarchy import linkage, dendrogram, cut_tree

class CustomParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def hierachichal(distancematrix):
    linkage_matrix = linkage(distancematrix, method='complete', metric='euclidean')
    print(linkage_matrix)
    cut = cut_tree(linkage_matrix, n_clusters = 3)
    print(cut)


def arguments():
    parser = CustomParser(description='Pairwise alignment (AVA).')
    parser.add_argument('--csv', type = str, required = True)
    #parser.add_argument('--outfile', type = str, required= True)
    
   
    return parser.parse_args()

def main():
    args = arguments()
    hierachichal(pd.read_csv(args.csv,index_col=str))


main()
