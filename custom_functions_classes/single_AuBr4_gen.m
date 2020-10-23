function AuBr4 = single_AuBr4_gen(AuBr4_coords,rel_pos,flat_CNT,steps,dopant_CNT_dist)
%single_AuBr4_gen Generates one AuBr4 molecule for band structure calcs
%	rel_pos = 0: centered on honeycomb, rel_pos = 1: centered on C atom;
%	unrolled_CNT = unrolled CNT object;
%	steps = array containing x_steps and z_steps.

d_CC = 1.421;	% bond length in graphite -- obtained from http://phycomp.technion.ac.il/~david/thesis/node3.html

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

x_steps = steps(1);
z_steps = steps(2);

if flat_CNT.chirality(1) == flat_CNT.chirality(2)	% armchair
	pos_x = d_CC*((3*x_steps-1)/2 - rel_pos);	% ref_atom x position
	pos_z = sqrt(3)*d_CC*z_steps/2;				% ref_atom z position
else	% zigzag
	pos_x = sqrt(3)*d_CC*(x_steps-1)/2;			% ref_atom x position
	pos_z = d_CC*((3*z_steps-1)/2 + rel_pos);	% ref_atom z position
end

pos_y = dopant_CNT_dist;
ref_atom_coords = [pos_x+flat_CNT.coord(1,1),(pos_y+flat_CNT.coord(1,2)),(pos_z+flat_CNT.coord(1,3))];

AuBr4 = molecule([gold,bromine_1,bromine_2,bromine_3,bromine_4],gold);
AuBr4.relative_pos(AuBr4_coords);
AuBr4.center_atom_duplication(ref_atom_coords);
AuBr4.rem_atoms_duplication();

end

