# Copyright 2020, Khai Yi Chin.
"""
Script to compute the RMS difference between coordinates of two .xyz files.
"""
from custom_classes_and_functions import compute_rms_difference, extract_xyz_file, XyzFileEntry

import argparse
import os

parser = argparse.ArgumentParser(description='Compute the RMS difference between coordinates of two .xyz files.')
parser.add_argument('file_1', type=str, help='The path to the 1st .xyz file.')
parser.add_argument('file_2', type=str, help='The path to the 2nd .xyz file.')
# parser.add_argument('--axis', type=str, help='The specific axis to evaluate.')

args = parser.parse_args()

# Populate arguments
file_1 = args.file_1
file_2 = args.file_2

# Read in both .xyz files
full_num_atoms, full_list_coord_obj_1 =  extract_xyz_file(os.path.abspath(file_1))
_, full_list_coord_obj_2 =  extract_xyz_file(os.path.abspath(file_2))

# Restructure object variables into lists of same axis
x_1 = [obj.x for obj in full_list_coord_obj_1]
x_2 = [obj.x for obj in full_list_coord_obj_2]

y_1 = [obj.y for obj in full_list_coord_obj_1]
y_2 = [obj.y for obj in full_list_coord_obj_2]

z_1 = [obj.z for obj in full_list_coord_obj_1]
z_2 = [obj.z for obj in full_list_coord_obj_2]

# Compute RMSD
print('\nComputing the RMSD values...')
rmsd_x = compute_rms_difference(x_1, x_2)
rmsd_y = compute_rms_difference(y_1, y_2)
rmsd_z = compute_rms_difference(z_1, z_2)

# Display output
print('\nSample size (total number of atoms per axis): ' + str(full_num_atoms))
print('\nRMSD x: ' + '{:.6f}'.format(rmsd_x) + ' angstroms')
print('\nRMSD y: ' + '{:.6f}'.format(rmsd_y) + ' angstroms')
print('\nRMSD z: ' + '{:.6f}'.format(rmsd_z) + ' angstroms')