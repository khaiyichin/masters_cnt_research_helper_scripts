function K = K_gen(K_coords,ref_atom_coords)
%K_gen Generates K atoms
%	K_coords = coordinates of the K atom, which is at the origin.
%	ref_atom_coords = array of atom coordinates, i.e. coordinates of K 
%	atoms.

potassium = dopant;
potassium.symbol = 'K';

K = molecule(potassium,potassium);
K.relative_pos(K_coords);
K.center_atom_duplication(ref_atom_coords);

end