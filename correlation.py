#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 13:05:11 2024

@author: fb445
"""
import sys
import numpy as np
from sklearn.cluster import SpectralClustering
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore", category=UserWarning, module="sklearn.manifold")
    
    
def readfile(filename):
    try:
        f = open(filename)
        out = f.readlines()
        f.close()
    except IOError:
        print('File %s does not exist!' % (filename))
        sys.exit()
    return out

def list2string(input_list,spacing=False,brackets=False, MCTDH = False):
    string = ''
    if brackets:
        string+='['
    if not spacing:
        for i in input_list:
            if MCTDH:
                string += 'Q' + str(i)
            else:
                string += str(i)
    else: 
        index = 0
        for i in input_list:
            if index == len(input_list)-1:
                if MCTDH:
                    string += 'Q' + str(i)
                else:
                    string += str(i)
            else:
                if MCTDH:
                    string += 'Q' + str(i)+" "
                else:
                    string += str(i)+" "
            index += 1
    if brackets:        
        string += "]"  
    
    return string          

def extract_submatrix(matrix, row_indices, col_indices):
    return matrix[np.ix_(row_indices, col_indices)]

def find(element, matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if matrix[i][j] == element:
                return (i, j)

# def correlationscore(Headmatrix, cluster):
    
#     score = 0
#     for mode_index, mode in enumerate(Headmatrix):
#         for compared_mode_index, compared_mode in enumerate(mode):
#             if mode_index in cluster and compared_mode_index in cluster:
#                 score += compared_mode
#             elif (mode_index not in cluster and compared_mode_index in cluster) or (mode_index in cluster and compared_mode_index not in cluster):
#                 score -= compared_mode
#     print("Score: " + str(score))        
#     return


class Vertex():

    def __init__(self, matrix, label):
        self.cluster_matrix = matrix
        self.labels = label
        self.clusters = []
        return

    def Clusterfunction(self,number_of_clusters, Min_number_of_clusters, depth = 0):
        # to check if the matrix should be divided into 2 clusters or more depending on the depth, after that spectral clustering is performed
        if number_of_clusters - depth >= Min_number_of_clusters:
            sc = SpectralClustering(number_of_clusters - depth, affinity='precomputed', n_init=100, random_state=69420)
            sc.fit(self.cluster_matrix)
        else: 
            sc = SpectralClustering(Min_number_of_clusters, affinity='precomputed', n_init=100, random_state=69420)
            sc.fit(self.cluster_matrix)
            
        All_clusters = []
        # spectral clustering yields the labeling and this is now sorted to find the matching matrix elements or normal mode labels
        for cluster_index in range(number_of_clusters - depth if number_of_clusters - depth >= Min_number_of_clusters else Min_number_of_clusters):
            cluster = []
            for item_index, item in enumerate(sc.labels_):
                if item == cluster_index:
                    cluster.append(self.labels[item_index])
            print("Cluster Nr. " + str(cluster_index+1) + ": " + list2string([i for i in cluster], spacing=True))
            All_clusters.append(cluster)
        
    
        All_cluster_matrices = []
        # now cluster matrix is divided according to the labels obtained before
        for cluster_index, cluster in enumerate(All_clusters):
            
            indicies = self.findlabels(cluster)

            sorted_matrix = extract_submatrix(self.cluster_matrix, indicies, indicies)
            # we build a new vertex and add it to self.clusters, after that clusterfunction is called again (recursion) for this subclusters
            new_vertex = Vertex(sorted_matrix, cluster)
            All_cluster_matrices.append(new_vertex)
            if len(new_vertex.labels) >= 3:
                new_depth = depth + 1
                new_vertex.Clusterfunction(number_of_clusters, Min_number_of_clusters, new_depth)
           
        self.clusters = All_cluster_matrices
        
        return  
        
    def findlabels(self, items):
            
        indicies = []
        for item in items:
            indicies.append(self.labels.index(item))
                
        return indicies
    
    # print the correlation matrix at the given vertex
    def plot(self):
            
        plt.imshow(self.cluster_matrix, cmap='viridis', interpolation='nearest')
        plt.colorbar()
        plt.xticks([])
        plt.yticks([])
        plt.gca().invert_yaxis()
        #plt.title("Heatmap of Correlation ")
       
        plt.savefig('correlation_matrix.eps', format='eps', bbox_inches='tight')
        plt.show()   
        
        return
    
    # return structure of the tree    
    def returnTree(self, depth=0, layers=None):

        if layers is None:
            layers = []

        if len(layers) <= depth:
            layers.append([])

        layers[depth].append(self.labels)
        for cluster in self.clusters:
            cluster.plotTree(depth + 1, layers)

        return layers
    
    # more simple plot of the tree not in the mctdh format
    def plotTree(self, depth=0, tree_info=None):
        
        if tree_info is None:
            tree_info = []
        
        indent = "    " * depth
        cluster_label = list2string(self.labels, spacing=True, brackets=True)
        tree_info.append(f"{indent}Cluster at depth {depth}: {cluster_label}")
        print(f"{indent}Cluster at depth {depth}: {cluster_label}")
        
        
        for subcluster in self.clusters:
            subcluster.plotTree(depth + 1, tree_info)
            
        return tree_info
    
    # returns the max depth of the tree
    def findMAXdepth(self, depth = 0, maxDepth = None):
        
        if maxDepth == None:
            maxDepth = 0
        
        if maxDepth < depth:
            maxDepth = depth
            
        for subcluster in self.clusters:
            maxDepth = subcluster.findMAXdepth(depth + 1, maxDepth)    
            
        return maxDepth
        
    
    def plotMCTDH(self, depth=0, tree_info=None, Nbr_of_SPF=str(2), initial = True):
        
        # Header of the mctdh input
        if tree_info is None:
            print("\n---Exemplary MCTDH Input---\n")
            print("ML-basis-section\n")
            tree_info = ["ML-basis-section\n"]
        # indentation depends on the depth of the tree
        indent = "\t" * depth
        # if there are not more subclusters -> plot the modes
        if not self.clusters:
            cluster_label = list2string(self.labels, spacing=True, brackets=True,MCTDH=True)
            tree_info.append(f"{indent}{depth}> {cluster_label}")
            print(f"{indent}{depth}> {cluster_label}")
        # If there are still subclusters we print the number of spf with the indentation. special case for the first layer as there is also the electronic branch
        else:
            if initial:
                cluster_label = " ".join(Nbr_of_SPF for i in range(0,len(self.clusters)+1))
                tree_info.append(f"{indent}{depth}> {cluster_label}")
                print(f"{indent}{depth}> {cluster_label}")
                
                depth2 = depth + 1
                indent2 = "\t" * depth2
                tree_info.append(f"{indent2}{depth+1}> [el]")
                print(f"{indent2}{depth+1}> [el]")
                
            else:
                cluster_label = " ".join(Nbr_of_SPF for i in range(0,len(self.clusters)))
                tree_info.append(f"{indent}{depth}> {cluster_label}")
                print(f"{indent}{depth}> {cluster_label}")
        # recursive call of this function for all subclusters (lower vertices)
        for subcluster in self.clusters:
            subcluster.plotMCTDH(depth + 1, tree_info, Nbr_of_SPF, initial = False)
        # final note at the end
        if depth == 0:
            tree_info.append("\n end-mlbasis-section")
            print("\nend-mlbasis-section")    
        
        return tree_info    

        
def main():
    
    # if len(sys.argv) != 4:
    #     print('Usage:\n python correlation.py <Matrix> <Max Number of Clusters> <Min Number of Clusters>\n')
    #     sys.exit()
    
    # lines = readfile(sys.argv[1])
    # Max_number_of_clusters = int(sys.argv[2])
    # Min_number_of_clusters = int(sys.argv[3])
  
    # Read the input file and how many clusters the matrix should be divided into e.g. Max = 4 and Min = 2 will lead to a seperation of the original matrix
    # into 4 clusters in the first layer, followed by 3 in the second layer and 2 in all remaining.
    lines = readfile('corr.dat')
    Max_number_of_clusters = 2
    Min_number_of_clusters = 2
    
    # Entries below the threshold will be set to 0
    threshold = 0
    
    print(f"--- Starting to Cluster for {Max_number_of_clusters} - {Min_number_of_clusters} Clusters ---")
    
    original_matrix = []
    
    for line in lines:
        original_matrix.append(line.split())
    
    original_matrix = np.array(original_matrix)
    original_matrix = original_matrix.astype(float)
    
    original_matrix = np.abs(original_matrix)
    
    # original_matrix = (np.ones(original_matrix.shape) + original_matrix)/2
    
    # W = C - 1
    weight_matrix = original_matrix - np.diag(np.ones(original_matrix.shape[0]))
    
    
    weight_matrix = np.abs(weight_matrix)
    weight_matrix[weight_matrix < threshold] = 0
    
    # Labels of the used correlation matrix, needs to have the same length as dimension of correlation matrix
    # indices = [7 ,8 ,9 ,11 ,12 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21 ,22 ,23 ,24 ,25 ,26 ,27 ,28 ,31 ,35 ,36 ,38 ,40 ,41 ,46 ,66 ,71 ,72 ,190 ,191]
    indices = [ 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 40, 41, 43, 45, 48, 49, 50, 54, 55, 57, 61, 66, 67, 68, 69, 71, 72, 77, 86, 96, 113, 114, 119, 133, 134, 135, 136, 142, 143, 144, 145, 148, 151, 161, 181, 188, 190, 191]

    # indices =  [n for n in range(1, len(weight_matrix)+1)]
    
    # build first vertex
    Headnote = Vertex(weight_matrix, indices)

    # Start building the Tree
    Headnote.Clusterfunction(Max_number_of_clusters, Min_number_of_clusters)
    
    print("Max Depth:" + str(Headnote.findMAXdepth()))
    
    Headnote.plot()
    
    # plot the tree in the mctdh format, also other ways of plotting are available (see above)
    tree_info = Headnote.plotMCTDH()

    return
    
if __name__ == '__main__':
         main()    

