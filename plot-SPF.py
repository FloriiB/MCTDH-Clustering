# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 15:36:43 2025

@author: flori
"""

import matplotlib.pyplot as plt
import numpy as np


def readfile(filename):
    try:
        f = open(filename)
        out = f.readlines()
        f.close()
    except IOError:
        print('File %s does not exist!' % (filename))
        exit()
    return out

def readMCTDHInp(file):

    lines = readfile(file)
    
    in_section = False
    mlbasis_section = []
    
    mctdh_start = 0
    mctdh_end = 0
    
    
    for line in lines:
        # expect input format like this: tfinal =2000.0 tout =50.0 with not space after the =
        if 'tfinal' in line:
            parts = line.split()
            for i, part in enumerate(parts):
                if part == 'tfinal':  
                    tfinal_value = float(parts[i + 1].replace('=', ''))
                    tout_value = float(parts[i + 3].replace('=', ''))
                break
         
    for line_index, line in enumerate(lines):
        if line.strip() == "ML-basis-section":
            mctdh_start = line_index
            in_section = True 
        elif line.strip() == "end-mlbasis-section":
            in_section = False
            mctdh_end = line_index
            mlbasis_section.append(line)
            break
        if in_section:
            mlbasis_section.append(line)    

    return lines, tfinal_value, tout_value, mlbasis_section, mctdh_start, mctdh_end

def findmaxDepth(input_lines):
    
    maxDepth = 0
    depth = 0
    
    for line in input_lines:
        if '>' in line:
            depth = int(line.split('>')[0])
        if depth > maxDepth:
            maxDepth = depth
        
    
    return maxDepth

def process_text(input_lines, maxDepth):
    
    numbers = [[] for _ in range(maxDepth + 1)]

    for line in input_lines:
        if '>' in line and 'el' not in line and 'Q' not in line:
            depth = int(line.split('>')[0])
            numbers[depth].extend(map(int, line.split('>')[1].strip().split()))

    return numbers

def read_SPF(mctdh_output, numbers, timesteps):
    
    size = 0
    for line in numbers:
        size += len(line)
    
    SPF_over_time = []
    searching_entry = True
    SPF_continue_reading = False
    
    index = 0
    depth = 0
    
    for line in mctdh_output:
        if ' fs,' in line and searching_entry:
            SPF_over_time.append([])
            searching_entry = False
            continue
        if 'layer' in line and not searching_entry:
            if depth != int(line.split()[-1]):
                depth = int(line.split()[-1])
                index = 0
                
            continue
        if line.startswith(' m') and not searching_entry:
            index += 1
            SPFs = line.split()[1:]
            length = len(SPFs)
            
            if length == numbers[depth][index-1]:
                SPF_over_time[-1].append(float(SPFs[-1]))
            elif length < numbers[depth][index-1]:
                SPF_continue_reading = True
                continue
        if SPF_continue_reading and not searching_entry:
            SPFs = line.split()
            length += len(SPFs)
            
            if length == numbers[depth][index-1]:
                SPF_over_time[-1].append(float(SPFs[-1]))
            else:
                continue
            SPF_continue_reading = False
            continue
        if 'Mode expectation values' in line:
            searching_entry = True
            continue
        
        
    return SPF_over_time
    
        

def main():
    
    inputfile  = 'mctdh.inp-threshold0'
    outputfile = 'output-threshold0'
    
    
    lines, tfinal_value, tout_value, mlbasis_section, mctdh_start, mctdh_end = readMCTDHInp(inputfile)
    numbers_by_layer = process_text(mlbasis_section, findmaxDepth(mlbasis_section))
    
    SPF_over_time = read_SPF(readfile(outputfile), numbers_by_layer, int(tfinal_value/tout_value))
    
    time = list(range(0,int(tfinal_value)+int(tout_value), int(tout_value)))
    
    SPF_over_time_transposed = list(map(list, zip(*SPF_over_time)))
    
    
    Averages = [np.mean(SPF[100:]) for SPF in SPF_over_time_transposed]
    
    print(np.mean(Averages))
    
    std = np.std(Averages)
    
    print(std)
    
    # for item_idx, item in enumerate(Averages):
    #     print("SPF Nr " + str(item_idx) + ": "+ str(item))
    # print("Standard Deviation: " + str(std))
    
    for Y in SPF_over_time_transposed:
        plt.plot(time[100:], Y[100:])
    plt.xlabel('time in fs')
    plt.ylabel('Coefficient of last SPF') 
    plt.yscale('log')
    
    
    #plt.savefig('coefficients.eps',format='eps',bbox_inches='tight')
    
    return
    
if __name__ == '__main__':
         main()   
