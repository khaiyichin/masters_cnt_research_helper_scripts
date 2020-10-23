function KAuBr4_coords = KAuBr4_atoms_avg_relaxed_coord(file)
%KAuBr4_atoms_avg_relaxed_coord Generates the averaged coordinates for all
%the KAuBr4 dopant atoms
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

K_index = find(species_num == 1);
K_x = x(K_index);
K_y = y(K_index);
K_z = z(K_index);

row_1_K_index = find(abs(K_y-row_1_pos)<3);
K_x = K_x(row_1_K_index);
K_y = K_y(row_1_K_index);
K_z = K_z(row_1_K_index);

Br_index = find(species_num == 3);
Br_x = x(Br_index);
Br_y = y(Br_index);
Br_z = z(Br_index);

row_1_Br_index = find(abs(Br_y-row_1_pos)<1.5);
Br_x = Br_x(row_1_Br_index);
Br_y = Br_y(row_1_Br_index);
Br_z = Br_z(row_1_Br_index);

total_dopants = length(row_1_ref_atom_index);
Br_atoms = 4;							% excluding the ref atom

for i = 1:total_dopants
	current_ind = (i-1)*Br_atoms+1;
	rel_K_x(i) = K_x(i) - single_row_ref_atom_coords(i,1);
	rel_K_y(i) = K_y(i) - single_row_ref_atom_coords(i,2);
	rel_K_z(i) = K_z(i) - single_row_ref_atom_coords(i,3);	
	
	rel_Br_x(current_ind:i*Br_atoms) = Br_x(current_ind:i*Br_atoms) - single_row_ref_atom_coords(i,1);
	rel_Br_y(current_ind:i*Br_atoms) = Br_y(current_ind:i*Br_atoms) - single_row_ref_atom_coords(i,2);
	rel_Br_z(current_ind:i*Br_atoms) = Br_z(current_ind:i*Br_atoms) - single_row_ref_atom_coords(i,3);
end

rel_K_x = rel_K_x';
rel_K_y = rel_K_y';
rel_K_z = rel_K_z';

rel_Br_x = rel_Br_x';
rel_Br_y = rel_Br_y';
rel_Br_z = rel_Br_z';

rel_Br_x = rel_Br_x';
rel_Br_y = rel_Br_y';
rel_Br_z = rel_Br_z';

% For KAuBr4 or AuBr4 dopants
SW = [];	% Br_1
NW = [];	% Br_2
NE = [];	% Br_3
SE = [];	% Br_4

for j = 1:total_dopants
	for k = 1:4
		current_ind = (j-1)*Br_atoms + k;

		% top
		if sign(rel_Br_x(current_ind)) == 1		% east
			if sign(rel_Br_z(current_ind)) == 1	% north
				NE = [NE;[rel_Br_x(current_ind),rel_Br_y(current_ind),rel_Br_z(current_ind)]];
			else								% south
				SE = [SE;[rel_Br_x(current_ind),rel_Br_y(current_ind),rel_Br_z(current_ind)]];
			end

		else									% west
			if sign(rel_Br_z(current_ind)) == 1	% north
				NW = [NW;[rel_Br_x(current_ind),rel_Br_y(current_ind),rel_Br_z(current_ind)]];
			else								% south
				SW = [SW;[rel_Br_x(current_ind),rel_Br_y(current_ind),rel_Br_z(current_ind)]];
			end
		end
	end
end

avg_SW = mean(SW);
avg_NW = mean(NW);
avg_NE = mean(NE);
avg_SE = mean(SE);
avg_K = [mean(rel_K_x),mean(rel_K_y),mean(rel_K_z)];

KAuBr4_coords = [zeros(1,3);avg_K;avg_SW;avg_NW;avg_NE;avg_SE];

end

