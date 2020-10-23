function K_coords = K_atoms_avg_relaxed_coord(file)
%K_atoms_avg_relaxed_coord Generates the averaged coordinates for all
%the K dopant atoms
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

dopant_index = find(species_num == 1);
x = x(dopant_index);
y = y(dopant_index);
z = z(dopant_index);

row_1_pos = y(1);
row_1_dop_atom_index = find(abs(y-row_1_pos)<1);
x = x(row_1_dop_atom_index);
y = y(row_1_dop_atom_index);
z = z(row_1_dop_atom_index);

K_coords = [mean(x);mean(y);mean(z)];

end

