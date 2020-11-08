from math import sqrt

"""Carbon-carbon bond length in angstroms"""
a = 1.42

"""Zigzag CNT periodic image spacing"""
zig_spacing = 1 * a

"""Armchair CNT periodic image spacing"""
arm_spacing = sqrt(3)/2 * a

"""Default x and y lattice spacing in angstroms"""
lattice_spacing = 20

"""Order of species in the job script's ChemicalSpeciesLabel block"""
species_dict = {'K': 1, 'Au': 2, 'Br': 3, 'C': 4}