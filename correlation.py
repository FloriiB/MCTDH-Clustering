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
        
        if number_of_clusters - depth >= Min_number_of_clusters:
            sc = SpectralClustering(number_of_clusters - depth, affinity='precomputed', n_init=100, random_state=69420)
            sc.fit(self.cluster_matrix)
        else: 
            sc = SpectralClustering(Min_number_of_clusters, affinity='precomputed', n_init=100, random_state=69420)
            sc.fit(self.cluster_matrix)
            
        All_clusters = []

        for cluster_index in range(number_of_clusters - depth if number_of_clusters - depth >= Min_number_of_clusters else Min_number_of_clusters):
            cluster = []
            for item_index, item in enumerate(sc.labels_):
                if item == cluster_index:
                    cluster.append(self.labels[item_index])
            print("Cluster Nr. " + str(cluster_index+1) + ": " + list2string([i for i in cluster], spacing=True))
            All_clusters.append(cluster)
        
    
        All_cluster_matrices = []
        
        for cluster_index, cluster in enumerate(All_clusters):
            
            indicies = self.findlabels(cluster)

            sorted_matrix = extract_submatrix(self.cluster_matrix, indicies, indicies)
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
    
    
    def __print__(self):
            
        plt.imshow(self.cluster_matrix, cmap='viridis', interpolation='nearest')
        plt.colorbar()
        plt.xticks([])
        plt.yticks([])
        plt.title("Heatmap of Cluster")
        plt.show()
            
        return
        
    def returnTree(self, depth=0, layers=None):

        if layers is None:
            layers = []

        if len(layers) <= depth:
            layers.append([])

        layers[depth].append(self.labels)
        for cluster in self.clusters:
            cluster.plotTree(depth + 1, layers)

        return layers

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
    
    def findMAXdepth(self, depth = 0, maxDepth = None):
        
        if maxDepth == None:
            maxDepth = 0
        
        if maxDepth < depth:
            maxDepth = depth
            
        for subcluster in self.clusters:
            maxDepth = subcluster.findMAXdepth(depth + 1, maxDepth)    
            
        return maxDepth
        
    
    def plotMCTDH(self, depth=0, tree_info=None, Nbr_of_SPF=str(2), initial = True):
        
        if tree_info is None:
            print("\n---Exemplary MCTDH Input---\n")
            print("ML-basis-section\n")
            tree_info = ["ML-basis-section\n"]
        
        indent = "    " * depth
        
        if not self.clusters:
            cluster_label = list2string(self.labels, spacing=True, brackets=True,MCTDH=True)
            tree_info.append(f"{indent} {depth}> {cluster_label}")
            print(f"{indent} {depth}> {cluster_label}")
        else:
            if initial:
                cluster_label = " ".join(Nbr_of_SPF for i in range(0,len(self.clusters)+1))
                tree_info.append(f"{indent} {depth}> {cluster_label}")
                print(f"{indent} {depth}> {cluster_label}")
                
                depth2 = depth + 1
                indent2 = "    " * depth2
                tree_info.append(f"{indent2} {depth+1}> [el]")
                print(f"{indent2} {depth+1}> [el]")
                
            else:
                cluster_label = " ".join(Nbr_of_SPF for i in range(0,len(self.clusters)))
                tree_info.append(f"{indent} {depth}> {cluster_label}")
                print(f"{indent} {depth}> {cluster_label}")

        for subcluster in self.clusters:
            subcluster.plotMCTDH(depth + 1, tree_info, Nbr_of_SPF, initial = False)
        
        if depth == 0:
            tree_info.append("\n end-mlbasis-section")
            print("\nend-mlbasis-section")    
        
        return tree_info    

        
def main():
    
    if len(sys.argv) != 4:
        print('Usage:\n python correlation.py <Matrix> <Max Number of Clusters> <Min Number of Clusters>\n')
        sys.exit()
    
    lines = readfile(sys.argv[1])
    Max_number_of_clusters = int(sys.argv[2])
    Min_number_of_clusters = int(sys.argv[3])
  
    
    print(f"--- Starting to Cluster for {Max_number_of_clusters} - {Min_number_of_clusters} Clusters ---")
    
    original_matrix = []
    
    for line in lines:
        original_matrix.append(line.split())
    
    original_matrix = np.array(original_matrix)
    original_matrix = original_matrix.astype(float)
    
    weight_matrix = original_matrix - np.diag(np.ones(original_matrix.shape[0]))
    
    Headnote = Vertex(weight_matrix, [n for n in range(7, len(weight_matrix)+7)])
    
    Headnote.Clusterfunction(Max_number_of_clusters, Min_number_of_clusters)
    
    print("Max Depth:" + str(Headnote.findMAXdepth()))
    
    tree_info = Headnote.plotMCTDH()

    return
    
if __name__ == '__main__':
         main()    

