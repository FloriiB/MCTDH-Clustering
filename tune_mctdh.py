# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 15:36:43 2025

@author: flori
"""

import sys
import numpy as np
import re
import copy

def readfile(filename):
    try:
        f = open(filename)
        out = f.readlines()
        f.close()
    except IOError:
        print('File %s does not exist!' % (filename))
        sys.exit()
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

    return lines, tfinal_value, mlbasis_section, mctdh_start, mctdh_end

def findmaxDepth(input_lines):
    
    maxDepth = 0
    depth = 0
    
    for line in input_lines:
        if '>' in line:
            depth = int(line.split('>')[0])
        if depth > maxDepth:
            maxDepth = depth
        
    
    return maxDepth

# def process_text(input_lines, maxDepth):
#     # Initialize a list of lists to store numbers for each layer
#     numbers = [[] for _ in range(maxDepth + 1)]

#     for line in input_lines:
#         if '>' in line and 'el' not in line and 'Q' not in line:
#             depth = int(line.split('>')[0])  # Extract the layer depth
#             # Extract all numbers after '>' and store them
#             numbers[depth].extend(map(int, line.split('>')[1].strip().split()))

#     return numbers

def process_text(input_lines, maxDepth):
    
    numbers = [[] for _ in range(maxDepth + 1)]

    for line in input_lines:
        if '>' in line and 'el' not in line and 'Q' not in line:
            depth = int(line.split('>')[0])
            numbers[depth].extend(map(int, line.split('>')[1].strip().split()))

    return numbers

def generate_text(input_lines, numbers_by_layer):
   
    numbers_copy = [iter(layer) for layer in numbers_by_layer]

    def replace_numbers(line):
        if '>' in line and 'el' not in line and 'Q' not in line:
            depth = int(line.split('>')[0])
            old_numbers = line.split('>')[1].strip().split()
            new_numbers = [str(next(numbers_copy[depth])) for _ in old_numbers]
            
            indentation = "\t" * depth
            return f"{indentation}{depth}> {' '.join(new_numbers)}"
        return line

    
    return [replace_numbers(line) for line in input_lines]

def updateSPFs(mctdh_output, tfinal_value, numbers, threshold):
    
    formatted_tfinal_value = f"{tfinal_value:.2f}"
    
    searching_entry = True
    SPF_continue_reading = False
    depth = 0
    index = 0
    
    for line in mctdh_output:
        if str(formatted_tfinal_value) + ' fs' in line and searching_entry:
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
                if float(SPFs[-1]) > threshold:
                    numbers[depth][index-1] += 1
            elif length < numbers[depth][index-1]:
                SPF_continue_reading = True
                continue
        if SPF_continue_reading and not searching_entry:
            SPFs = line.split()
            length += len(SPFs)
            
            if length == numbers[depth][index-1]:
                if float(SPFs[-1]) > threshold:
                    numbers[depth][index-1] += 1
            else:
                continue
            
            SPF_continue_reading = False
            continue
    
    
    
    return numbers



def main():
    
    inputfile  = 'mctdh.inp'
    outputfile = 'output'
    
    threshold = 10
    
    
    lines, tfinal_value, mlbasis_section, mctdh_start, mctdh_end = readMCTDHInp(inputfile)
    
    
    
    numbers_by_layer = process_text(mlbasis_section, findmaxDepth(mlbasis_section))
    
    print('current SPFs:')
    for depth, numbers in enumerate(numbers_by_layer):
        print(f"Layer {depth}: {numbers}")
    print('\n')
    
    numbers_by_layer = updateSPFs(readfile(outputfile), tfinal_value, numbers_by_layer, threshold)
    
    print('changed SPFs:')
    for depth, numbers in enumerate(numbers_by_layer):
        print(f"Layer {depth}: {numbers}")
    print('\n')
    
    updated_text = generate_text(mlbasis_section, numbers_by_layer)
    
    lines[mctdh_start:mctdh_end] = updated_text
    
    
    for line in updated_text:
        print(line)
    
    
    
    
    
    
    return
    
if __name__ == '__main__':
         main()   