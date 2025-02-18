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
            return f"{indentation}{depth}> {' '.join(new_numbers)}\n"
        return line

    
    return [replace_numbers(line) for line in input_lines]

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

def update_SPFs(mctdh_output, numbers, threshold):
    
    
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
        if ' fs,' in line and searching_entry:
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
                if float(SPFs[-1]) > threshold and numbers[depth][index-1] != reference_numbers[depth][index-1]:
                    numbers[depth][index-1] += 1
                    changes_made = True
            elif length < numbers[depth][index-1]:
                SPF_continue_reading = True
                continue
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
        if 'Mode expectation values' in line:
            searching_entry = True
            continue
        
    return numbers, changes_made
    
def run():
    
    workdir = os.path.join(os.environ['SCRATCH_DIR'])
    
    os.popen('cp mctdh.inp ' + workdir+'/mctdh.inp')
    os.popen('cp mctdh.op ' + workdir+'/mctdh.op')
    
    time.sleep(3)
    
    string = 'mctdh86 -mnd mctdh'
    runerror = runProgram(string, workdir, 'output')
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
        time.sleep(3) # I need something better than this... We need to wait for the copy step to be finished before proceeding
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
    threshold = 5
    
    
    
    
    run_process = True
    
    while run_process:
    
        lines, tfinal_value, mlbasis_section, mctdh_start, mctdh_end = readMCTDHInp(inputfile)

        numbers_by_layer = process_text(mlbasis_section, findmaxDepth(mlbasis_section))
    
    
        print('current SPFs:')
        for depth, numbers in enumerate(numbers_by_layer):
            print(f"Layer {depth}: {numbers}")
        print('\n')
    
        run()
        if only_final_value:
            numbers_by_layer, run_process = update_final_SPFs(readfile(outputfile), tfinal_value, numbers_by_layer, threshold)
        else:
            numbers_by_layer, run_process = update_SPFs(readfile(outputfile), numbers_by_layer, threshold)
        
        print('changed SPFs:')
        for depth, numbers in enumerate(numbers_by_layer):
            print(f"Layer {depth}: {numbers}")
        print('\n')
    
    
        updated_text = generate_text(mlbasis_section, numbers_by_layer)
        
        lines[mctdh_start:mctdh_end+1] = updated_text
    
        writefile(inputfile, lines)
    
    
    return
    
if __name__ == '__main__':
         main()   
