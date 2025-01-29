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

def main():
    
    outputfile = 'output'
    operatorfile = 'mctdh.op'
    
    Emin = 1.6
    Emax = 2.8

    K = 30
    
    
    plt.rcParams['font.family'] = 'serif'  # Use 'serif' to allow Times New Roman
    plt.rcParams['font.serif'] = ['Times New Roman']  # Specify the font
    plt.rcParams.update({'font.size': 18})   
    
    
    
    NRG_S, NRG_T = read_operator(readfile(operatorfile))
    
    Population, Time = read_pop(readfile(outputfile))
    
    Sinlget_pop = [ entry[:len(NRG_S)]  for entry in Population]
    Triplet_pop = [ entry[len(NRG_S):]  for entry in Population]
    
    NRG_S.pop(0)
    
    # Prepare Singlet Pops
    Sinlget_pop = list(map(list, zip(*Sinlget_pop)))
    Sinlget_pop.pop(0)
    Sinlget_pop = list(map(list, zip(*Sinlget_pop)))
    Sinlget_pop = np.matrix(Sinlget_pop)
    
    Nbr_of_triplets = int(len(NRG_T)/3)
    for index, row in enumerate(Triplet_pop):
        help_triplet = [row[:Nbr_of_triplets], row[Nbr_of_triplets:2*Nbr_of_triplets], row[2*Nbr_of_triplets:]]
        help_triplet = list(map(list, zip(*help_triplet)))
        final = []
        for item in help_triplet:
            final += item
        Triplet_pop[index] = final
        
    
    Triplet_pop = np.matrix(Triplet_pop)
    
    Energies = np.linspace(Emin, Emax, 91)

    # Create empty arrays for both populations
    Z_full_S = np.zeros((len(Time), len(Energies)))
    Z_full_T = np.zeros((len(Time), len(Energies)))

    # Fill the Singlet population matrix
    used_energy_idx = []
    for idx, energy in enumerate(NRG_S):
        energy_idx = np.argmin(np.abs(Energies - energy))
        if energy_idx not in used_energy_idx:
            used_energy_idx.append(energy_idx)
            Z_full_S[:, energy_idx] = Sinlget_pop[:, idx].flatten()
        else:
            Z_full_S[:, energy_idx] = Z_full_S[:, energy_idx] + Sinlget_pop[:, idx].flatten()

    # Fill the Triplet population matrix
    used_energy_idx = []
    for idx, energy in enumerate(NRG_T):
        energy_idx = np.argmin(np.abs(Energies - energy))  
        if energy_idx not in used_energy_idx:
            used_energy_idx.append(energy_idx)
            Z_full_T[:, energy_idx] = Triplet_pop[:, idx].flatten()
        else: 
            Z_full_T[:, energy_idx] =  Z_full_T[:, energy_idx] + Triplet_pop[:, idx].flatten()
            
    # Transpose the matrices for plotting
    Z_S = Z_full_S.T  
    Z_T = Z_full_T.T  

    vmin = min(Z_S.min(), Z_T.min())
    vmax = max(Z_S.max(), Z_T.max())

    # Plot side by side
    fig, axes = plt.subplots(1, 2, figsize=(8, 6), gridspec_kw={'wspace': 0}, sharey=True)

    # Singlet plot
    im1 = axes[0].imshow(Z_S, aspect='auto', cmap='hot', origin='upper', vmin=vmin, vmax=vmax)
    axes[0].set_title('Singlet Population')
    
    # Triplet plot
    im2 = axes[1].imshow(Z_T, aspect='auto', cmap='hot', origin='upper', vmin=vmin, vmax=vmax)
    axes[1].set_title('Triplet Population')

    # X-axis ticks
    num_ticks_x = 5
    tick_indices_x = np.linspace(0, len(Time) - 1, num_ticks_x, dtype=int)
    tick_indices_x = tick_indices_x[:-1]
    for ax in axes:
        ax.set_xticks(tick_indices_x)
        ax.set_xticklabels([Time[i] for i in tick_indices_x])
        ax.set_xlabel('Time in fs')
        
    # Y-axis ticks
    rounded_labels = [f"{e:.2f}" for e in Energies[::K]]
    rounded_labels_help = [e for e in Energies[::K]]
    energy_tick_indices = [np.where(Energies == item)[0][0] for item in rounded_labels_help if item in Energies]
    
    axes[0].set_yticks(energy_tick_indices)
    axes[0].set_yticklabels(rounded_labels)
    axes[0].set_ylabel('Energy in eV')
    axes[0].set_ylim(ax.get_ylim()[::-1])
    
    axes[1].tick_params(axis='y', which='both', left=False, right=False)
    
    # Add single color bar
    cbar_ax = fig.add_axes([0.95, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im1, cax=cbar_ax, orientation='vertical')
    cbar.set_label('Population')

    #plt.savefig('singlet_triplet_comparison.png', format='png', bbox_inches='tight')
    plt.show()
    
    return
    
if __name__ == '__main__': 
    main()
    
    






















