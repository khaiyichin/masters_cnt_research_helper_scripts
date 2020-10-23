function AuBr4_coords = AuBr4_atoms_avg_relaxed_coord(file)
%AuBr4_atoms_avg_relaxed_coord Generates the averaged coordinates for all
%the AuBr4 dopant atoms
%   file = relaxed model .xyz file
[species,x,y,z] = importRelaxedxyz(file);

% K - 1; Au - 2; Br - 3; C - 4
species_num = zeros(length(species),1);
for num_of_atoms = 1:length(species)
	switch species{num_of_atoms}
		case 'K'
			species_num(num_of_atoms) = 1;
		case 'Au'
			species_num(num_of_atoms) = 2;
		case 'Br'
			species_num(num_of_atoms) = 3;
		case 'C'
			species_num(num_of_atoms) = 4;
	end
end

dopant_index = find(species_num ~= 4);
x = x(dopant_index);
y = y(dopant_index);
z = z(dopant_index);
species_num = species_num(dopant_index);

ref_atom_index = find(species_num == 2);		% Au atoms as reference
ref_atom_coords = [x(ref_atom_index) y(ref_atom_index) z(ref_atom_index)];
row_1_pos = max(ref_atom_coords(:,2));			% top row as row_1
row_1_ref_atom_index = find(abs(ref_atom_coords(:,2)-row_1_pos)<1);

single_row_ref_atom_coords = [ref_atom_coords(row_1_ref_atom_index,1)...
	ref_atom_coords(row_1_ref_atom_index,2) ref_atom_coords(row_1_ref_atom_index,3)];

non_ref_index = find(species_num ~= 2);
x = x(non_ref_index);
y = y(non_ref_index);
z = z(non_ref_index);

row_1_rem_atom_index = find(abs(y-row_1_pos)<2);
x = x(row_1_rem_atom_index);
y = y(row_1_rem_atom_index);
z = z(row_1_rem_atom_index);

total_dopants = length(row_1_ref_atom_index);
rem_dopant_atoms = 4;							% excluding the ref atom

for i = 1:total_dopants
	current_ind = (i-1)*rem_dopant_atoms+1;
	rel_x(current_ind:i*rem_dopant_atoms) = x(current_ind:i*rem_dopant_atoms) - single_row_ref_atom_coords(i,1);
	rel_y(current_ind:i*rem_dopant_atoms) = y(current_ind:i*rem_dopant_atoms) - single_row_ref_atom_coords(i,2);
	rel_z(current_ind:i*rem_dopant_atoms) = z(current_ind:i*rem_dopant_atoms) - single_row_ref_atom_coords(i,3);
end

rel_x = rel_x';
rel_y = rel_y';
rel_z = rel_z';

rel_x = rel_x';
rel_y = rel_y';
rel_z = rel_z';

% For KAuBr4 or AuBr4 dopants
SW = [];	% Br_1
NW = [];	% Br_2
NE = [];	% Br_3
SE = [];	% Br_4

for j = 1:total_dopants
	for k = 1:4
		current_ind = (j-1)*rem_dopant_atoms + k;

		% top
		if sign(rel_x(current_ind)) == 1		% east
			if sign(rel_z(current_ind)) == 1	% north
				NE = [NE;[rel_x(current_ind),rel_y(current_ind),rel_z(current_ind)]];
			else								% south
				SE = [SE;[rel_x(current_ind),rel_y(current_ind),rel_z(current_ind)]];
			end

		else									% west
			if sign(rel_z(current_ind)) == 1	% north
				NW = [NW;[rel_x(current_ind),rel_y(current_ind),rel_z(current_ind)]];
			else								% south
				SW = [SW;[rel_x(current_ind),rel_y(current_ind),rel_z(current_ind)]];
			end
		end
	end
end

avg_SW = mean(SW);
avg_NW = mean(NW);
avg_NE = mean(NE);
avg_SE = mean(SE);

AuBr4_coords = [zeros(1,3);avg_SW;avg_NW;avg_NE;avg_SE];

end

