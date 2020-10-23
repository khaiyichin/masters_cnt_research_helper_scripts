function K = single_K_gen(K_coords,rel_pos,flat_CNT,steps,dopant_CNT_dist)
%single_K_gen Generates one K molecule for band structure calcs
%	rel_pos = 0: centered on honeycomb, rel_pos = 1: centered on C atom;
%	unrolled_CNT = unrolled CNT object;
%	steps = array containing x_steps and z_steps.

d_CC = 1.421;	% bond length in graphite -- obtained from http://phycomp.technion.ac.il/~david/thesis/node3.html

potassium = dopant;
potassium.symbol = 'K';

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

K = molecule(potassium,potassium);
K.relative_pos(K_coords);
K.center_atom_duplication(ref_atom_coords);

end

