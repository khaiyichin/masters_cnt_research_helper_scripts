function AuBr4 = AuBr4_gen(AuBr4_coords,ref_atom_coords)
%AuBr4_gen Generates KAuBr4 molecules
%	AuBr4_coords = coordinates of the 5 atoms, with the 'center' atom
%	(reference atom), i.e. Au atom at the origin.
%	ref_atom_coords = array of 'center' atom coordinates, i.e. coordinates
%	of Au atoms.

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

AuBr4 = molecule([gold,bromine_1,bromine_2,bromine_3,bromine_4],gold);
AuBr4.relative_pos(AuBr4_coords);
AuBr4.center_atom_duplication(ref_atom_coords);
AuBr4.rem_atoms_duplication();

end