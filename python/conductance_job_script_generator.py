# Copyright 2020, Khai Yi Chin.
"""
Script to generate a conductance job script based on the .xyz file generated from a relaxation run.
"""
from custom_classes_and_functions import XyzFileEntry, ConductanceMain, ConductanceFullDevice, ConductanceElectrode

import argparse
import os

parser = argparse.ArgumentParser(description='Generate conductance job script based on .xyz file.')
parser.add_argument('filename', type=str, help='The path to the .xyz file.')
parser.add_argument('--output', type=str, help='The path to the generated script.')

args = parser.parse_args()

# Populate arguments
filename = args.filename
if not args.output:
    output = os.path.dirname(filename)
else: output = args.output

"""Prompt for SBATCH configuration"""
# Required inputs
while True:
    job_name = input('Job name: ')
    
    if job_name:
        break
    else:
        print('Job name is required!')

num_tot_dop = int(input('Number of dopant atoms in the model: '))
num_elec_atoms = int(input('Number of atoms in each electrode: '))
num_elec_dop = int(input('Number of dopant atoms in each electrode: ')) if num_tot_dop > 0 else 0
cnt_type = input('CNT type, zigzag or armchair? (zigzag): ') or 'zigzag'
tacc_type = input('TACC system type, ls5 or st2? (ls5): ') or 'ls5'

# Optional inputs
should_modify = input('Modify default values? (no): ') or 'no'

if should_modify == 'y' or should_modify == 'yes' or should_modify == 'Y':
    job_queue = input('Job queue name (normal): ') or 'normal'
    job_nodes = input('Number of nodes (2): ') or '2'
    job_ntasks = input('Number of tasks per node (24): ') or '24'
    job_runtime = input('Maximum allowed runtime in hh:mmm:ss (48:00:00): ') or '48:00:00'
    job_email = input('Notification email address (khaiyichin@utexas.edu): ') or 'khaiyichin@utexas.edu'
    job_notif_type = input('Notification types (all): ') or 'all'
    job_alloc = input('Allocation name, Multiscale-Simulatio or AIMD? (Multiscale-Simulatio): ') or 'Multiscale-Simulatio'
else:
    job_queue = 'normal'
    job_nodes = '2'
    job_ntasks = '24'
    job_runtime = '48:00:00'
    job_email = 'khaiyichin@utexas.edu'
    job_notif_type = 'all'
    job_alloc = 'Multiscale-Simulatio'

"""Script generation"""
# Generate the main conductance script
main = ConductanceMain(job_name, job_queue, job_nodes, job_ntasks, job_runtime, job_email, job_notif_type, job_alloc)
main.generate_script()

# Generate the full device job script
full_dev = ConductanceFullDevice(filename, job_name, num_elec_atoms, cnt_type, num_elec_dop, num_tot_dop, tacc_type)
full_dev.generate_script()

# Generate the electrodes
left_elec = ConductanceElectrode(filename, 'left', num_elec_atoms, cnt_type, num_elec_dop, num_tot_dop, tacc_type)
left_elec.generate_script()

# Generate the electrodes
right_elec = ConductanceElectrode(filename, 'right', num_elec_atoms, cnt_type, num_elec_dop, num_tot_dop, tacc_type)
right_elec.generate_script()