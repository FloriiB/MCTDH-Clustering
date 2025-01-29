# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 15:36:43 2025

@author: flori
"""

import sys
import os
import subprocess as sp
import time
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

def writefile(filename, content):
    try:
        f = open(filename, 'w')
        if isinstance(content, list):
            for line in content:
                f.write(line)
        elif isinstance(content, str):
            f.write(content)
        else:
            print('Content %s cannot be written to file!' % (content))
        f.close()
    except IOError:
        print('Could not write to file %s!' % (filename))
        sys.exit()
    return

#read MCTDH input and return all the lines, the final time values, the separate mlbasis section und the indices where to insert the last part
def readMCTDHInp(file):

    lines = readfile(file)
    
    in_section = False
    mlbasis_section = []
    
    mctdh_start = 0
    mctdh_end = 0
    
    
    for line in lines:
        # expect input format like this: tfinal =2000.0 tout =50.0 with not space after the =. Otherwise the program will crash, this should be improved
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

# this function reads to mlbasis section and returns all numbers of spf in a large list and returns it
def process_text(input_lines, maxDepth):
    
    numbers = [[] for _ in range(maxDepth + 1)]

    for line in input_lines:
        if '>' in line and 'el' not in line and 'Q' not in line:
            depth = int(line.split('>')[0])
            numbers[depth].extend(map(int, line.split('>')[1].strip().split()))

    return numbers

# takes the original input and a new set of numbers and generates a new mlbasis section
def generate_text(input_lines, numbers_by_layer):
   
    numbers_copy = [iter(layer) for layer in numbers_by_layer]

    def replace_numbers(line):
        if '>' in line and 'el' not in line and 'Q' not in line:
            depth = int(line.split('>')[0])
            old_numbers = line.split('>')[1].strip().split()
            new_numbers = [str(next(numbers_copy[depth])) for _ in old_numbers]
            
            indentation = "\t" * depth
            return f"{indentation}{depth}> {' '.join(new_numbers)}\n"
        return line

    
    return [replace_numbers(line) for line in input_lines]


# goes through the mctdh output and increases the number of spfs if the last spf coefficient is larger than threshold (detailed explanation see below)
def update_final_SPFs(mctdh_output, tfinal_value, numbers, threshold):
    
    formatted_tfinal_value = f"{tfinal_value:.2f}"
    
    changes_made = False
    
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
                    changes_made = True
            elif length < numbers[depth][index-1]:
                SPF_continue_reading = True
                continue
        if SPF_continue_reading and not searching_entry:
            SPFs = line.split()
            length += len(SPFs)
            
            if length == numbers[depth][index-1]:
                if float(SPFs[-1]) > threshold:
                    numbers[depth][index-1] += 1
                    changes_made = True
            else:
                continue
            
            SPF_continue_reading = False
            continue
    
    return numbers, changes_made

# same as above but we take into account all timestamps and not just the final one
def update_SPFs(mctdh_output, numbers, threshold):
    
    # we copy the current numbers of spf to compare them after the update. This is currently used to ensure we increase the amount of spfs only by +1
    # every iteration. One could rework this to take an increase by e.g. +2 or a decrease into account
    reference_numbers = copy.deepcopy(numbers)
    

    for i in range(len(reference_numbers)):
        for j in range(len(reference_numbers[i])):
            reference_numbers[i][j] += 1
    
    
    changes_made = False
    searching_entry = True
    SPF_continue_reading = False
    depth = 0
    index = 0
    
    for line in mctdh_output:
        # skip all lines until we find an entry marked by the timestamp
        if ' fs,' in line and searching_entry:
            searching_entry = False
            continue
        # take layer to find the depth
        if 'layer' in line and not searching_entry:
            if depth != int(line.split()[-1]):
                depth = int(line.split()[-1])
                index = 0
                
            continue
        
        # entries of the spfs coefficents start with " m" so we take this line
        if line.startswith(' m') and not searching_entry:
            index += 1
            SPFs = line.split()[1:]
            length = len(SPFs)
            # sometimes there is a linebreak if the amount of spf is too large. So we count the amount of spf in our currently line until
            # until the length (amount of spfs in the current line + previous lines) is equal to expected number predicted by the array "numbers"
            if length == numbers[depth][index-1]:
                # this if means "if the last spf coefficient is smaller than the threshold and if we havent done any changes to this number of spf so far"
                if float(SPFs[-1]) > threshold and numbers[depth][index-1] != reference_numbers[depth][index-1]:
                    numbers[depth][index-1] += 1
                    changes_made = True
            # if we didnt find the end yet we do a linebreak and keep reading untul length == numbers[depth][index-1]
            elif length < numbers[depth][index-1]:
                SPF_continue_reading = True
                continue
        # just like above, we keep reading until we find the end of the coefficients
        if SPF_continue_reading and not searching_entry:
            SPFs = line.split()
            length += len(SPFs)
            
            if length == numbers[depth][index-1]:
                if float(SPFs[-1]) > threshold and numbers[depth][index-1] != reference_numbers[depth][index-1]:
                    numbers[depth][index-1] += 1
                    changes_made = True
            else:
                continue
            SPF_continue_reading = False
            continue
        # end of the timestamp, we turn searching_entry back on to skip lines until we find another entry
        if 'Mode expectation values' in line:
            searching_entry = True
            continue
        
    return numbers, changes_made
    
def run():
    # workdir needs to be set in the slurm file, one could also define one right here
    workdir = os.path.join(os.environ['SCRATCH_DIR'])
    
    os.popen('cp mctdh.inp ' + workdir+'/mctdh.inp')
    os.popen('cp mctdh.op ' + workdir+'/mctdh.op')
    
    time.sleep(3) # I need something better than this... We need to wait for the copy step to be finished before proceeding. Similar issue later on
    # call mctdh
    string = 'mctdh86 -mnd mctdh'
    runerror = runProgram(string, workdir, 'output')
    
    # program will exit in the case of a crash of the mctdh calculation
    if runerror.returncode != 0:
        print('Seems like MCTDH is in severe trouble...')
        sys.exit()

    return

def runProgram(string, workdir, outfile, errfile=''):
    
    prevdir = os.getcwd()
    os.chdir(workdir)

    try:
        with open(os.path.join(workdir, outfile), 'w') as stdoutfile:
            if errfile:
                with open(os.path.join(workdir, errfile), 'w') as stderrfile:
                    runerror = sp.run(string, shell=True, stdout=stdoutfile, stderr=stderrfile)
            else:
                runerror = sp.run(string, shell=True, stdout=stdoutfile, stderr=sp.STDOUT)
    except OSError as e:
        print('Call had some serious problems:', e)
        sys.exit()
    finally:
        os.popen('cp mctdh/output ' + prevdir)
        
        print('MCTDH Calculation finished!')
        time.sleep(3)
        os.popen('ls')
        print('cleaning up scratch directory: ' + os.getcwd() + '\n')
        os.popen('rm -r mctdh')
        time.sleep(3)
        os.popen('rm *')
        time.sleep(3)
        os.chdir(prevdir)
        
    return runerror

def main():
    
    inputfile  = 'mctdh.inp'
    outputfile = 'output'
    only_final_value = False
    threshold = 10
    
    
    
    
    run_process = True
    
    while run_process:
        
        # read input
        lines, tfinal_value, mlbasis_section, mctdh_start, mctdh_end = readMCTDHInp(inputfile)

        # extract number of spf from mlbasis section
        numbers_by_layer = process_text(mlbasis_section, findmaxDepth(mlbasis_section))
    
        
        print('current SPFs:')
        for depth, numbers in enumerate(numbers_by_layer):
            print(f"Layer {depth}: {numbers}")
        print('\n')
    
        # run program and update number of spf based on mctdh output
        run()
        if only_final_value:
            numbers_by_layer, run_process = update_final_SPFs(readfile(outputfile), tfinal_value, numbers_by_layer, threshold)
        else:
            numbers_by_layer, run_process = update_SPFs(readfile(outputfile), numbers_by_layer, threshold)
        
        print('changed SPFs:')
        for depth, numbers in enumerate(numbers_by_layer):
            print(f"Layer {depth}: {numbers}")
        print('\n')
    
        # generate new mlbasis section with new numbers of spf
        updated_text = generate_text(mlbasis_section, numbers_by_layer)
        
        # insert new mlbasis section
        lines[mctdh_start:mctdh_end+1] = updated_text
    
        # write mctdh input
        writefile(inputfile, lines)
    
    
    return
    
if __name__ == '__main__':
         main()   
