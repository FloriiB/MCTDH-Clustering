# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 15:36:43 2025

@author: flori
"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def readfile(filename):
    try:
        f = open(filename)
        out = f.readlines()
        f.close()
    except IOError:
        print('File %s does not exist!' % (filename))
        exit()
    return out


def read_pop(lines):
    
    reading = False
    Population = []
    Time = []
    
    for line in lines:
        if "population" in line:
            reading = True
            
            entries = line.split()[2:]
            entries = [float(item) for item in entries]
            Population.append(entries)
            
            continue
        
        if "Q" in line:
            reading = False
            continue
        if reading:
            entries = line.split()
            entries = [float(item) for item in entries]
            Population[-1] = Population[-1] + entries
        if "Time" in line:
            Time.append(float(line.split()[2]))
        
    return Population, Time
    
def read_operator(lines):
    
    NRG_S = []
    NRG_T = []
    
    
    for line in lines:
        
        if line.startswith('eS'):
            try:
                NRG_S.append(float(line.split('=')[1].split(',')[0]))
            except:
                a=1
        if line.startswith('eT'):
            try:
                NRG_T.append(float(line.split('=')[1].split(',')[0]))
                NRG_T.append(float(line.split('=')[1].split(',')[0]))
                NRG_T.append(float(line.split('=')[1].split(',')[0]))
            except:
                a=1
    return  NRG_S, NRG_T    

def remove_duplicates(lst):
    seen = set()
    result = []
    for item in lst:
        if item not in seen:
            seen.add(item)
            result.append(item)
    return result

def calculatePop(Sinlget_pop, Triplet_pop):
    
    MLCT_pop = []
    MC_pop = []
    Groundstate_pop = []
    
    for line_index, Singlet_line in enumerate(Sinlget_pop):
        
        
        MLCT_current_pop = 0
        MC_current_pop = 0
        
        for Singlet_item_index, Singlet_item in enumerate(Singlet_line):
            if Singlet_item_index in Singlet_Labels_MC:
                MC_current_pop += Singlet_item
            elif Singlet_item_index in Singlet_Labels_MLCT:
                MLCT_current_pop += Singlet_item
        
        MLCT_pop.append(MLCT_current_pop)
        MC_pop.append(MC_current_pop)
        Groundstate_pop.append(Singlet_line[0])
        
        
          
    for line_index, Triplet_line in enumerate(Triplet_pop):
            
       Triplet_ms0 = Triplet_line[0::3]
       Triplet_ms1 = Triplet_line[1::3]
       Triplet_ms1m = Triplet_line[2::3]
       
       MLCT_current_pop = 0
       MC_current_pop = 0
       
       for item_index, item in enumerate(Triplet_ms0):
           if item_index+1 in Triplet_Labels_MC:
               MC_current_pop += Triplet_ms0[item_index]
               MC_current_pop += Triplet_ms1[item_index]
               MC_current_pop += Triplet_ms1m[item_index]
               
           elif item_index+1 in Triplet_Labels_MLCT:
               MLCT_current_pop += Triplet_ms0[item_index]
               MLCT_current_pop += Triplet_ms1[item_index]
               MLCT_current_pop += Triplet_ms1m[item_index]
        
       MLCT_pop[line_index] += MLCT_current_pop
       MC_pop[line_index] += MLCT_current_pop
    
    
    return MLCT_pop, MC_pop

def readInputs(outputfile, operatorfile):
    
    NRG_S, NRG_T = read_operator(readfile(operatorfile))
    
    Population, Time = read_pop(readfile(outputfile))
    
    Sinlget_pop = [ entry[:len(NRG_S)]  for entry in Population]
    Triplet_pop = [ entry[len(NRG_S):]  for entry in Population]

    
    Nbr_of_triplets = int(len(NRG_T)/3)
    for index, row in enumerate(Triplet_pop):
        help_triplet = [row[:Nbr_of_triplets], row[Nbr_of_triplets:2*Nbr_of_triplets], row[2*Nbr_of_triplets:]]
        help_triplet = list(map(list, zip(*help_triplet)))
        final = []
        for item in help_triplet:
            final += item
        Triplet_pop[index] = final
        
    
    return Time, Sinlget_pop, Triplet_pop

import csv
def lists_to_csv(list1, list2, filename="output.csv"):
    # Ensure both lists have the same length
    max_length = max(len(list1), len(list2))
    list1.extend([None] * (max_length - len(list1)))
    list2.extend([None] * (max_length - len(list2)))
    
    # Write to CSV file
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Column 1", "Column 2"])  # Header row
        for item1, item2 in zip(list1, list2):
            writer.writerow([item1, item2])
    
    print(f"CSV file '{filename}' has been created successfully.")

def main():
    
    #outputfiles = ['output-threshold0-3', 'output-threshold1-2', 'output-threshold2-2', 'output-threshold3', 'output-threshold3-5', 'output']
    outputfiles = ['output-20', 'output-10', 'output-5', 'reference']
    operatorfile = 'mctdh.op'
    
    global Singlet_Labels_MC, Singlet_Labels_MLCT, Triplet_Labels_MC, Triplet_Labels_MLCT
    Singlet_Labels_MC = [7,8,9]
    Singlet_Labels_MLCT = [1,2,3,4,5,6]
    Triplet_Labels_MC = [1,2,3,10,11,12]
    Triplet_Labels_MLCT = [4,5,6,7,8,9]
    
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman']
    plt.rcParams.update({'font.size': 18})   
    
    custom_labels = ["20", "10", "5", "Reference"]
    
    cmap = plt.get_cmap("prism")

    num_files = len(outputfiles)
    colors = [cmap(i / num_files) for i in range(num_files)]  # Generate evenly spaced colors
    
    for i, outputfile in enumerate(outputfiles):
    
        Time, Sinlget_pop, Triplet_pop = readInputs(outputfile, operatorfile)
    
        MLCT_pop, MC_pop = calculatePop(Sinlget_pop, Triplet_pop)
    
        plt.plot(Time, MLCT_pop, label = custom_labels[i], color = colors[i])
        #plt.plot(Time, MC_pop, color='b')
    #plt.plot(Time, Groundstate_pop)
    plt.legend(loc="upper right")
    plt.xlabel('Time [fs]')
    plt.ylabel('MLCT Population')
    plt.ylim(-0.01, 0.28)    
    plt.savefig('Threshold1-pop.eps', format='eps', bbox_inches='tight')
    
    #plt.xlim(0, 1000)
    return
    
if __name__ == '__main__': 
    main()
    
    






















