function KAuBr4 = KAuBr4_gen(KAuBr4_coords,ref_atom_coords)
%KAuBr4_gen Generates KAuBr4 molecules
%	KAuBr4_coords = coordinates of the 6 atoms, with the 'center' atom
%	(reference atom), i.e. Au atom at the origin.
%	ref_atom_coords = array of 'center' atom coordinates, i.e. coordinates
%	of Au atoms.

potassium = dopant;
potassium.symbol = 'K';

gold = dopant;
gold.symbol = 'Au';

bromine_1 = dopant;
bromine_1.symbol = 'Br';

bromine_2 = dopant;
bromine_2.symbol = 'Br';

bromine_3 = dopant;
bromine_3.symbol = 'Br';

bromine_4 = dopant;
bromine_4.symbol = 'Br';

KAuBr4 = molecule([gold,potassium,bromine_1,bromine_2,bromine_3,bromine_4],gold);
KAuBr4.relative_pos(KAuBr4_coords);
KAuBr4.center_atom_duplication(ref_atom_coords);
KAuBr4.rem_atoms_duplication();

end