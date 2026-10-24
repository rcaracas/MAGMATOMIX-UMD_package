#!/usr/bin/env python3
# coding: utf-8
### Ana Anzulović

import os, re, subprocess, sys
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages
import csv
import matplotlib.cm as cm

def find_files(current_dir,r):
    umd_files = []
    bonding_files = []
    for filename in os.listdir(current_dir):
        if filename.endswith('.outcar.umd.dat'):
            umd_files.append(filename)
        elif filename.endswith(f'{r}.popul.dat'):
            bonding_files.append(filename)
    matching_pairs = []
    for umd_file in umd_files:
        common_prefix = extract_common_prefix([umd_file])
        matching_bonding_files = [file for file in bonding_files if file.startswith(common_prefix)]
        for matching_bonding_file in matching_bonding_files:
            matching_pairs.append((umd_file, matching_bonding_file))
    n = len(matching_pairs)
    return matching_pairs, n

def extract_common_prefix(file_list):
    common_prefix = os.path.commonprefix(file_list)
    return re.sub(r'\.umd\.dat$|\.bonding\.dat$', '', common_prefix)

def read_umd(file_path):
    density_list, pressure_list, temperature_list, acell, elements = [], [], [], [], []
    simulation_time = 0
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith('Density'):
                density_list.append(float(line.split()[1]))
            elif line.startswith('Pressure'):
                pressure_list.append(float(line.split()[1]))
            elif line.startswith('Temperature'):
                temperature_list.append(float(line.split()[1]))
            elif line.startswith('time'):
                simulation_time = max(simulation_time, float(line.split()[1]))
            elif line.startswith('acell'):
                acell.append(float(line.split()[1]))
            elif line.startswith('elements'):
                elements.extend(line.split()[1:])
    return {'Density': np.array(density_list), 'Pressure': np.array(pressure_list), 'Temperature': np.array(temperature_list), 'time': simulation_time, 'acell': acell, 'elements': elements}

def read_populfile(file_path):
    clusters_data = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines[4:]:
            data = line.strip().split()
            cluster_name = data[0]
            cluster_time = float(data[3])
            atomic_types = re.findall(r'([A-Za-z]+)_\d+', cluster_name)
            atomic_numbers = [int(num) for num in re.findall(r'(\d+)', cluster_name)]
            clusters_data.append({'cluster_name': cluster_name, 'cluster_time': cluster_time, 'atomic_types': atomic_types, 'atomic_numbers': atomic_numbers})
    return clusters_data

def find_and_filter_clusters(clusters_data, file_name, pressure_mean, atomic_type_flags):
    """
    Filters clusters based on the number of flags provided and writes them to lifetimes.txt.
    - 1 flag: Finds clusters made of ONLY that element.
    - 2+ flags: Finds clusters containing ALL of the specified elements.
    """
    num_flags = len(atomic_type_flags)
    
    with open("lifetimes.txt", "a") as output_file:
        for cluster in clusters_data:
            cluster_name = cluster['cluster_name']
            atomic_types_in_cluster = set(cluster['atomic_types'])
            
            # Mode 1: Single-element search (e.g., find H_2, H_3...)
            if num_flags == 1:
                flag = atomic_type_flags[0]
                if len(atomic_types_in_cluster) == 1 and flag in atomic_types_in_cluster:
                    cluster_info = [cluster_name, str(cluster['cluster_time']), str(cluster['atomic_types']), str(cluster['atomic_numbers'])]
                    output_file.write(f"{file_name} {pressure_mean} {' '.join(cluster_info)}\n")
            
            # Mode 2: Multi-flag search (e.g., find clusters with both O and H)
            elif num_flags >= 2:
                # Check if all required flags are present in the cluster
                if all(flag in atomic_types_in_cluster for flag in atomic_type_flags):
                    cluster_info = [cluster_name, str(cluster['cluster_time']), str(cluster['atomic_types']), str(cluster['atomic_numbers'])]
                    output_file.write(f"{file_name} {pressure_mean} {' '.join(cluster_info)}\n")

def grep_pattern(FileName, Pattern, SkipSteps,simulation_time): 
    data, average, stdev, variance = [], 0, 0, 0
    anchor=Pattern.split()[0]
    try:
        patterns=subprocess.check_output(['grep',Pattern,FileName]).decode("utf-8")
    except subprocess.CalledProcessError:
        return [], 0, 0
    greps=patterns.split('\n')
    for isteps in range(SkipSteps,len(greps)):
        if not greps[isteps]: continue
        elems=greps[isteps].split()
        for ii, elem in enumerate(elems):
            if elem == anchor and ii + 1 < len(elems):
                try: data.append(float(elems[ii+1])); break
                except ValueError: continue
    if not data: return [], 0, 0
    average = sum(data)/len(data)
    variance = sum([(x-average)**2 for x in data]) / len(data)
    stdev = np.sqrt(variance)
    return data, average, stdev

def get_cluster_composition(cluster_name):
    parts = re.findall(r'([A-Za-z]+)_(\d+)', cluster_name)
    return {atom: int(num) for atom, num in parts}

def get_canonical_key(composition_dict):
    """Creates a unique, sorted key from a composition dict, e.g., (('H', 1), ('O', 1))"""
    return tuple(sorted(composition_dict.items()))

def save_plotting_data_to_csv(data_dict, output_file):
    with open(output_file, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Cluster", "Pressure Mean (GPa)", "Lifetime (fs)"])
        for pressure, clusters in data_dict.items():
            for cluster, lifetimes in clusters.items():
                for lifetime in lifetimes:
                    writer.writerow([cluster, pressure, lifetime])

def main():
    if '-h' in sys.argv or '--help' in sys.argv:
        help_text = f"""
    DESCRIPTION:
    This script analyzes cluster lifetimes from molecular dynamics simulation data.
    It processes UMD and POPUL files, calculates lifetimes for specific cluster
    compositions, and generates a grid of square plots saved to a PDF.

    FILE REQUIREMENTS:
    The script must be run in a directory containing pairs of simulation files:
        - `r<r_value>.umd.dat`
        - `r<r_value>.popul.dat`

    USAGE:
    lifetimes.py <r_value> <flag1> [<flag2> ...]

    ARGUMENTS:
    <r_value>         An integer corresponding to the r value (e.g., 1 or 0).
    <flag1>           The chemical symbol of the element to search for (e.g., 'H', 'O').
    [<flag2>]         (Optional) The chemical symbol of the second element.
    """
    print(help_text)
    sys.exit(0)
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <r_value> <flag1> [<flag2>]")
        print("\nMODES:")
        print(f"  1. Single-Element Search: {sys.argv[0]} 1 H")
        print(f"  2. Multi-Flag Search: {sys.argv[0]} 1 O H")
        print(f"     -> Finds clusters containing BOTH O and H, e.g., O_1H_1, C_1O_2H_1...")
        sys.exit(1)
    
    r = int(sys.argv[1])
    
    # Get all remaining arguments from the command line as a list of flags
    atomic_type_flags = sys.argv[2:]
    
    num_flags = len(atomic_type_flags)
    print(f"Searching for lifetimes of clusters containing {', '.join(sys.argv[2:])}.")

    current_dir = os.getcwd()
    matching_pairs, _ = find_files(current_dir,r)
    if os.path.exists("lifetimes.txt"): os.remove("lifetimes.txt")
        
    for umd_file, bonding_file in matching_pairs:
        umd_data = read_umd(os.path.join(current_dir, umd_file))
        clusters_data = read_populfile(os.path.join(current_dir, bonding_file))
        _, pressure_mean, _ = grep_pattern(os.path.join(current_dir, umd_file), 'Pressure', 0, umd_data['time'])
        find_and_filter_clusters(clusters_data, bonding_file, pressure_mean, atomic_type_flags)

    second_entry_dict = defaultdict(lambda: defaultdict(list))
    with open("lifetimes.txt", "r") as file:
        for line in file:
            parts = line.split()
            if len(parts) < 4: continue
            second_entry, cluster, fourth_entry = round(float(parts[1]), 2), parts[2], float(parts[3])
            second_entry_dict[second_entry][cluster].append(fourth_entry)

    if not second_entry_dict:
        print("No lifetime data found for the given flags. Exiting.")
        return

    csv_filename = f"plotting_data_r{r}_{'_'.join(atomic_type_flags)}.csv"
    save_plotting_data_to_csv(second_entry_dict, csv_filename)
    print(f"Raw plotting data saved to {csv_filename}")
    
    sorted_second_entries = sorted(second_entry_dict.keys())
    all_clusters = set(c for data in second_entry_dict.values() for c in data.keys())

    unique_row_keys = set()
    if num_flags == 1:
        # Mode 1: Key is just the count of the single atom
        flag1 = atomic_type_flags[0]
        for cluster in all_clusters:
            composition = get_cluster_composition(cluster)
            unique_row_keys.add((composition.get(flag1, 0),))
    elif num_flags >= 2:
        # Mode 2: Key is the tuple of counts of the first two flags
        flag1, flag2 = atomic_type_flags[0], atomic_type_flags[1]
        for cluster in all_clusters:
            composition = get_cluster_composition(cluster)
            unique_row_keys.add((composition.get(flag1, 0), composition.get(flag2, 0)))
    
    sorted_row_keys = sorted(list(unique_row_keys))
    if not sorted_row_keys:
        print("No clusters matching the criteria were found. Cannot generate plot.")
        return

    colormap = plt.colormaps['viridis'].resampled(len(sorted_row_keys))
    with PdfPages(f"lifetimes_r{r}_{'_'.join(atomic_type_flags)}.pdf") as pdf:
        columns, rows = len(sorted_second_entries), len(sorted_row_keys)
        # 1. Define the desired final appearance of ONE subplot and the spacing
        subplot_width = 5.0
        subplot_height = 5.0
        subplot_spacing = {
        'left':   0.21,
        'right':  0.90,
        'bottom': 0.16,
        'top':    0.92,
        'wspace': 0.2,
        'hspace': 0.4
    }

        # 2. Calculate the total figure size needed to accommodate this
        total_horizontal_plot_space = subplot_width * (columns + subplot_spacing['wspace'] * (columns - 1))
        figure_width = total_horizontal_plot_space / (subplot_spacing['right'] - subplot_spacing['left'])

        total_vertical_plot_space = subplot_height * (rows + subplot_spacing['hspace'] * (rows - 1))
        figure_height = total_vertical_plot_space / (subplot_spacing['top'] - subplot_spacing['bottom'])

        # 3. Create the figure and axes with the calculated size
        fig, axes = plt.subplots(
            rows, 
            columns, 
            figsize=(figure_width, figure_height), # Use the calculated dimensions
            squeeze=False, 
            sharey=False
        )

        # 4. Apply the spacing adjustments
        fig.subplots_adjust(**subplot_spacing)   
        for col_idx, second_entry in enumerate(sorted_second_entries):
            cluster_data = second_entry_dict.get(second_entry, {})
            for row_idx, key_tuple in enumerate(sorted_row_keys):
                ax = axes[row_idx, col_idx]
                matching_clusters = []
                for cluster, fourth_entries in cluster_data.items():
                    composition = get_cluster_composition(cluster)
                    # Check if the cluster's composition matches the row key for the current mode
                    current_key = ()
                    if num_flags == 1:
                        current_key = (composition.get(atomic_type_flags[0], 0),)
                    elif num_flags >= 2:
                        current_key = (composition.get(atomic_type_flags[0], 0), composition.get(atomic_type_flags[1], 0))
                    
                    if current_key == key_tuple:
                        matching_clusters.append((cluster, fourth_entries))

                if matching_clusters:
                    offset, all_lifetimes_in_cell = 0, []
                    matching_clusters.sort(key=lambda x: x[0])
                    
                    for cluster, fourth_entries in matching_clusters:
                        sorted_lifetimes = sorted(fourth_entries, reverse=True)
                        all_lifetimes_in_cell.extend(sorted_lifetimes)
                        bar_pos = np.arange(offset, offset + len(sorted_lifetimes)) + 1
                        ax.bar(bar_pos, sorted_lifetimes, color=colormap(row_idx), alpha=0.7, label=cluster)
                        offset += len(sorted_lifetimes)
                    
                    avg_lifetime = sum(all_lifetimes_in_cell) / len(all_lifetimes_in_cell)
                    ax.text(0.95, 0.95, f'Avg: {avg_lifetime:.0f} fs', transform=ax.transAxes, fontsize=40, ha='right', va='top')
                    ax.tick_params(axis='both', which='major', labelsize=30)
                    
                    if row_idx == 0: ax.set_title(f'{second_entry} GPa', fontsize=30)
                    if col_idx == 0: ax.set_ylabel('Lifetime (fs)', fontsize=30)
                    if row_idx == rows - 1: ax.set_xlabel('Occurrences', fontsize=30)

                    if col_idx == columns - 1:
                        subscript_map = str.maketrans('0123456789', '₀₁₂₃₄₅₆₇₈₉')
                        label_text = ""
                        if num_flags == 1:
                            label_text = f"{atomic_type_flags[0]}{str(key_tuple[0]).translate(subscript_map)}"
                        elif num_flags >= 2:
                            flag1, flag2 = atomic_type_flags[0], atomic_type_flags[1]
                            label_text = f"{flag1}{str(key_tuple[0]).translate(subscript_map)}{flag2}{str(key_tuple[1]).translate(subscript_map)}"
                        
                        ax.text(1.05, 0.5, label_text, transform=ax.transAxes, fontsize=30, color='k', va='center', rotation=90)
                    
                    if len(matching_clusters) > 1: ax.legend(title="Isomers", fontsize=10)
                else:
                    ax.axis('off')
        
        pdf.savefig(fig, transparent=True)
        plt.show()
    

if __name__ == "__main__":
    main()