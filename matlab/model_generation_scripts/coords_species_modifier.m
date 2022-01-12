%% Write to text file with imported coordinates and modified species number
% Only converts .xyz files into .txt files
clear

cd /home/khaiyichin/masters_cnt_research/

folder_to_read = 'temp/';
file_to_read = [folder_to_read,'80SWNT_K_1.xyz'];

folder_to_write = '/home/khaiyichin/masters_cnt_research/temp/test';
file_to_write = [folder_to_write,'80SWNT_K_1.xyz'];
if ~exist(folder_to_write, 'dir')
  mkdir(folder_to_write);
end

%% Check whether reading .xyz file of relaxed coordinates (CHECK BEFORE RUNNING)
relaxed_coords = 3;		% 1 - unrelaxed (coordinates for relaxation calcs); 2 - unrelaxed (coordinates for conduction calcs); 3 - relaxed; 4 - Band structure calcs
big_box = 0;			% 1 = big computational unit cell; 0 = translational asymmetric unit cell in z direction
chirality = [8,0];

% For relaxed_coords ~= 1 or 4, fill the parameters below:
num_of_CNT = 1;
electrode_num_atoms = 64;	% C only without dopants in ONE electrode: for SWNT (5,5) each unit cell has 20 atoms; for (8,0) each unit cell has 32 atoms
electrode_dopant_atoms = 2;	% total dopant atoms in ONE electrode

%% Convert coordinates and species
% Accounting for semi-infinity and periodicity of electrodes
d_CC  = 1.421;							% bond length in graphite -- obtained from http://phycomp.technion.ac.il/~david/thesis/node3.html
if chirality(1) == chirality(2)			% armchair config
	some_space = d_CC*sqrt(3)/2;
end

if	chirality(2) == 0					% zig-zag config
	some_space = d_CC;
end

if relaxed_coords == 3
	[species,x,y,z] = importRelaxedxyz(file_to_read);		% SIESTA generates xyz files with space delimiters, which is different from the xyz file generator we used (tab delimited)
else
	[species,x,y,z] = importxyz(file_to_read);
end

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

switch relaxed_coords 
	case 1			% pre-relaxed coord transformation
		
		% Sorting C atoms first
		C_index = find(species_num == 4);
		x_C = x(C_index);
		y_C = y(C_index);
		z_C = z(C_index);
		len_x = max(x_C) - min(x_C);			% circumferential length
		len_y = max(y_C) - min(y_C);
		len_z = max(z_C) - min(z_C);			% axial length
		
		CNT_center_x = (max(x_C) + min(x_C))/2;
		CNT_center_y = (max(y_C) + min(y_C))/2;
		CNT_center_z = (max(z_C) + min(z_C))/2;

		cell_size_x = len_x + 40;			% for an additional 20 on each side, + 40;
		cell_size_y = len_y + 40;
		
		if big_box == 0
			cell_size_z = len_z + some_space;
			z = z - z(1);
		else
			cell_size_z = len_z + 40;
			z = z + (cell_size_z/2 - CNT_center_z);
		end

		x = x + (cell_size_x/2 - CNT_center_x);
		y = y + (cell_size_y/2 - CNT_center_y);
		
		dopant_index = find(species_num ~= 4);		% find only the dopant molecule atoms

		combined_mat = [x(dopant_index),y(dopant_index),z(dopant_index),species_num(dopant_index)];
		sorted_mat = sortrows(combined_mat,3);		% sort rows based on the values in the 3rd column (z-coords)

		x(dopant_index) = sorted_mat(:,1);
		y(dopant_index) = sorted_mat(:,2);
		z(dopant_index) = sorted_mat(:,3);
		species_num(dopant_index) = sorted_mat(:,4);

		fileID = fopen(file_to_write,'w');
		fprintf(fileID,'%s\t%s\t%s\t%s\t90.\t90.\t90.\n\n',['from ',file_to_read,'; num of atoms = ',num2str(num_of_atoms),'; lattice param: '],...
			num2str(cell_size_x),num2str(cell_size_y),num2str(cell_size_z));
		fprintf(fileID,'%s\t\t%s\t\t%s\t\t%s\n','x','y','z','chemical species number');
		fprintf(fileID,'%4.5f\t%4.6f\t%4.6f\t %2i\r\n',[x,y,z,species_num]');
		fclose(fileID);
		
	case 2			% unrelaxed coord transformation (for unrelaxed conduction calculation)
		dopant_index = find(species_num ~= 4);		% find only the dopant molecule atoms

		combined_mat = [x(dopant_index),y(dopant_index),z(dopant_index),species_num(dopant_index)];
		sorted_mat = sortrows(combined_mat,3);		% sort rows based on the values in the 3rd column (z-coords)

		x(dopant_index) = sorted_mat(:,1);
		y(dopant_index) = sorted_mat(:,2);
		z(dopant_index) = sorted_mat(:,3);
		species_num(dopant_index) = sorted_mat(:,4);

		z_carbon = z(1:length(z)-length(dopant_index));	% only the C atoms
		atoms_per_layer = chirality(1) + chirality(2);
		atoms_per_layer = num_of_CNT*atoms_per_layer;
		left = [];
		right = [];

		for i = 1:atoms_per_layer			% obtaining the two end layers
			[value_1,index_1] = min(z_carbon);
			z_carbon(index_1) = NaN;
			[value_2,index_2] = max(z_carbon);
			z_carbon(index_2) = NaN;

			left = [left;value_1,index_1];
			right = [right;value_2,index_2];		
		end

		len_x = max(x) - min(x);
		len_y = max(y) - min(y);
		len_of_CNT = max(z) - min(z);

		ave_pos_right = mean(right(:,1));		% averaging the end layer positions (so each end layer share the same z-coord)
		ave_pos_left = mean(left(:,1));

		cell_size_x = len_x + 40;
		cell_size_y = len_y + 40;	
		cell_size_z = (ave_pos_right - ave_pos_left) + some_space;	% the space is to ensure the electrodes are semi-infinite and periodic, i.e. the xy plane of the rightside of the box consists of (first layer) atoms from the repeated CNT (periodic BCs)

		CNT_center_x = (max(x) + min(x))/2;
		CNT_center_y = (max(y) + min(y))/2;

		x = x + (cell_size_x/2 - CNT_center_x);
		y = y + (cell_size_y/2 - CNT_center_y);
		z = z - ave_pos_left;
		z(left(:,2)) = 0;						% making the leftmost atoms having the same z-coords
		z(right(:,2)) = ave_pos_right - ave_pos_left;	% making the rightmost atoms having the same z-coords

		% Arranging into sections of coordinates
		last_right_ind = max(right(:,2));		% last C atom index (used for the right electrode)

		% Left electrode
		left_electrode_x = x(1:electrode_num_atoms);
		left_electrode_y = y(1:electrode_num_atoms);
		left_electrode_z = z(1:electrode_num_atoms);
		left_electrode_species = species_num(1:electrode_num_atoms);
		last_left_C_atom_coord = min(left_electrode_z(end-sum(chirality)+1:end));
		if electrode_dopant_atoms ~= 0
			left_electrode_x = [left_electrode_x;x(last_right_ind+1:last_right_ind+electrode_dopant_atoms)];
			left_electrode_y = [left_electrode_y;y(last_right_ind+1:last_right_ind+electrode_dopant_atoms)];
			left_electrode_z = [left_electrode_z;z(last_right_ind+1:last_right_ind+electrode_dopant_atoms)];
			left_electrode_species = [left_electrode_species;species_num(last_right_ind+1:last_right_ind+electrode_dopant_atoms)];
		end

		% Right electrode
		right_electrode_x = x(last_right_ind-electrode_num_atoms+1:last_right_ind);
		right_electrode_y = y(last_right_ind-electrode_num_atoms+1:last_right_ind);
		right_electrode_z = z(last_right_ind-electrode_num_atoms+1:last_right_ind);
		right_electrode_species = species_num(last_right_ind-electrode_num_atoms+1:last_right_ind);
		last_right_C_atom_coord = min(right_electrode_z(end-sum(chirality)+1:end));
		first_right_C_atom_coord = min(right_electrode_z(1:sum(chirality)));
		if electrode_dopant_atoms ~= 0
			right_electrode_x = [right_electrode_x;x(end-electrode_dopant_atoms+1:end)];
			right_electrode_y = [right_electrode_y;y(end-electrode_dopant_atoms+1:end)];
			right_electrode_z = [right_electrode_z;z(end-electrode_dopant_atoms+1:end)];
			right_electrode_species = [right_electrode_species;species_num(end-electrode_dopant_atoms+1:end)];
		end

		% Device
		device_x = x(electrode_num_atoms+1:last_right_ind-electrode_num_atoms);
		device_y = y(electrode_num_atoms+1:last_right_ind-electrode_num_atoms);
		device_z = z(electrode_num_atoms+1:last_right_ind-electrode_num_atoms);
		device_species = species_num(electrode_num_atoms+1:last_right_ind-electrode_num_atoms);
		if electrode_dopant_atoms ~= 0
			device_x = [device_x;x(last_right_ind+electrode_dopant_atoms+1:end-electrode_dopant_atoms)];
			device_y = [device_y;y(last_right_ind+electrode_dopant_atoms+1:end-electrode_dopant_atoms)];
			device_z = [device_z;z(last_right_ind+electrode_dopant_atoms+1:end-electrode_dopant_atoms)];
			device_species = [device_species;species_num(last_right_ind+electrode_dopant_atoms+1:end-electrode_dopant_atoms)];
		else
			device_x = [device_x;x(last_right_ind+1:end)];
			device_y = [device_y;y(last_right_ind+1:end)];
			device_z = [device_z;z(last_right_ind+1:end)];
			device_species = [device_species;species_num(last_right_ind+1:end)];
		end

		full_model_x = [left_electrode_x;device_x;right_electrode_x];
		full_model_y = [left_electrode_y;device_y;right_electrode_y];
		full_model_z = [left_electrode_z;device_z;right_electrode_z];
		full_model_species = [left_electrode_species;device_species;right_electrode_species];
		
		fileID = fopen(file_to_write,'w');	
		fprintf(fileID,'%s\n\n',['from ',file_to_read]);
		fprintf(fileID,'%s\t\t%s\t\t%s\t\t%s\n','x','y','z','chemical species number');
		fprintf(fileID,'\n%s-atom left electrode with lattice param: %s\t%s\t%s\t90.\t90.\t90.\n',...
			num2str(electrode_num_atoms+electrode_dopant_atoms),num2str(cell_size_x),num2str(cell_size_y),num2str(last_left_C_atom_coord+some_space));
		fprintf(fileID,'%4.5f\t%4.5f\t%4.5f\t%2i\r\n',[left_electrode_x,left_electrode_y,left_electrode_z,left_electrode_species]');
		fprintf(fileID,'\n%s-atom right electrode with lattice param: %s\t%s\t%s\t90.\t90.\t90.\n',...
			num2str(electrode_num_atoms+electrode_dopant_atoms),num2str(cell_size_x),num2str(cell_size_y),num2str(last_right_C_atom_coord-first_right_C_atom_coord+some_space));	
		fprintf(fileID,'%4.5f\t%4.5f\t%4.5f\t%2i\r\n',[right_electrode_x,right_electrode_y,right_electrode_z,right_electrode_species]');
		fprintf(fileID,'\n%s-atom full model (left+device+right) with lattice param: %s\t%s\t%s\t90.\t90.\t90.\n',num2str(num_of_atoms),num2str(cell_size_x),...
			num2str(cell_size_y),num2str(cell_size_z));
		fprintf(fileID,'%4.5f\t%4.5f\t%4.5f\t%2i\r\n',[full_model_x,full_model_y,full_model_z,full_model_species]');
		fclose(fileID);

	case 3			% post-relaxed coord transformation (only changing the z-coordinates)
	dopant_index = find(species_num ~= 4);		% find only the dopant molecule atoms
	
	combined_mat = [x(dopant_index),y(dopant_index),z(dopant_index),species_num(dopant_index)];
	sorted_mat = sortrows(combined_mat,3);		% sort rows based on the values in the 3rd column (z-coords)
	
	x(dopant_index) = sorted_mat(:,1);
	y(dopant_index) = sorted_mat(:,2);
	z(dopant_index) = sorted_mat(:,3);
	species_num(dopant_index) = sorted_mat(:,4);
	
% 	z_coord_array = z;
	z_carbon = z(1:length(z)-length(dopant_index));	% only the C atoms
	atoms_per_layer = chirality(1) + chirality(2);
	atoms_per_layer = num_of_CNT*atoms_per_layer;
	left = [];
	right = [];
	
	for i = 1:atoms_per_layer			% obtaining the two end layers
		[value_1,index_1] = min(z_carbon);
		z_carbon(index_1) = NaN;
		[value_2,index_2] = max(z_carbon);
		z_carbon(index_2) = NaN;
		
		left = [left;value_1,index_1];
		right = [right;value_2,index_2];		
	end
		
	len_x = max(x) - min(x);
	len_y = max(y) - min(y);
    max(y)
    min(y)
	len_of_CNT = max(z) - min(z);

	ave_pos_right = mean(right(:,1));		% averaging the end layer positions (so each end layer share the same z-coord)
	ave_pos_left = mean(left(:,1));
	
	cell_size_x = len_x + 40;
	cell_size_y = len_y + 40;	
	cell_size_z = (ave_pos_right - ave_pos_left) + some_space;	% the space is to ensure the electrodes are semi-infinite and periodic, i.e. the xy plane of the rightside of the box consists of (first layer) atoms from the repeated CNT (periodic BCs)

	z = z - ave_pos_left;
	z(left(:,2)) = 0;						% making the leftmost atoms having the same z-coords
	z(right(:,2)) = ave_pos_right - ave_pos_left;	% making the rightmost atoms having the same z-coords
	
	% Arranging into sections of coordinates
	last_right_ind = max(right(:,2));		% last C atom index (used for the right electrode)
	
	% Left electrode
	left_electrode_x = x(1:electrode_num_atoms);
	left_electrode_y = y(1:electrode_num_atoms);
	left_electrode_z = z(1:electrode_num_atoms);
	left_electrode_species = species_num(1:electrode_num_atoms);
	last_left_C_atom_coord = min(left_electrode_z(end-sum(chirality)+1:end));
% 	last_left_C_atom_coord = left_electrode_z(end);
	if electrode_dopant_atoms ~= 0
		left_electrode_x = [left_electrode_x;x(last_right_ind+1:last_right_ind+electrode_dopant_atoms)];
		left_electrode_y = [left_electrode_y;y(last_right_ind+1:last_right_ind+electrode_dopant_atoms)];
		left_electrode_z = [left_electrode_z;z(last_right_ind+1:last_right_ind+electrode_dopant_atoms)];
		left_electrode_species = [left_electrode_species;species_num(last_right_ind+1:last_right_ind+electrode_dopant_atoms)];
	end
	
	% Right electrode
	right_electrode_x = x(last_right_ind-electrode_num_atoms+1:last_right_ind);
	right_electrode_y = y(last_right_ind-electrode_num_atoms+1:last_right_ind);
	right_electrode_z = z(last_right_ind-electrode_num_atoms+1:last_right_ind);
	right_electrode_species = species_num(last_right_ind-electrode_num_atoms+1:last_right_ind);
	last_right_C_atom_coord = min(right_electrode_z(end-sum(chirality)+1:end));
	first_right_C_atom_coord = min(right_electrode_z(1:sum(chirality)));
	if electrode_dopant_atoms ~= 0
		right_electrode_x = [right_electrode_x;x(end-electrode_dopant_atoms+1:end)];
		right_electrode_y = [right_electrode_y;y(end-electrode_dopant_atoms+1:end)];
		right_electrode_z = [right_electrode_z;z(end-electrode_dopant_atoms+1:end)];
		right_electrode_species = [right_electrode_species;species_num(end-electrode_dopant_atoms+1:end)];
	end
	
	% Device
	device_x = x(electrode_num_atoms+1:last_right_ind-electrode_num_atoms);
	device_y = y(electrode_num_atoms+1:last_right_ind-electrode_num_atoms);
	device_z = z(electrode_num_atoms+1:last_right_ind-electrode_num_atoms);
	device_species = species_num(electrode_num_atoms+1:last_right_ind-electrode_num_atoms);
	if electrode_dopant_atoms ~= 0
		device_x = [device_x;x(last_right_ind+electrode_dopant_atoms+1:end-electrode_dopant_atoms)];
		device_y = [device_y;y(last_right_ind+electrode_dopant_atoms+1:end-electrode_dopant_atoms)];
		device_z = [device_z;z(last_right_ind+electrode_dopant_atoms+1:end-electrode_dopant_atoms)];
		device_species = [device_species;species_num(last_right_ind+electrode_dopant_atoms+1:end-electrode_dopant_atoms)];
	else
		device_x = [device_x;x(last_right_ind+1:end)];
		device_y = [device_y;y(last_right_ind+1:end)];
		device_z = [device_z;z(last_right_ind+1:end)];
		device_species = [device_species;species_num(last_right_ind+1:end)];
	end
	
	full_model_x = [left_electrode_x;device_x;right_electrode_x];
	full_model_y = [left_electrode_y;device_y;right_electrode_y];
	full_model_z = [left_electrode_z;device_z;right_electrode_z];
	full_model_species = [left_electrode_species;device_species;right_electrode_species];
	
	fileID = fopen(file_to_write,'w');	
	fprintf(fileID,'%s\n\n',['from ',file_to_read]);
	fprintf(fileID,'%s\t\t%s\t\t%s\t\t%s\n','x','y','z','chemical species number');
	fprintf(fileID,'\n%s-atom left electrode with lattice param: %s\t%s\t%s\t90.\t90.\t90.\n',...
		num2str(electrode_num_atoms+electrode_dopant_atoms),num2str(cell_size_x),num2str(cell_size_y),num2str(last_left_C_atom_coord+some_space));
	fprintf(fileID,'%4.5f\t%4.5f\t%4.5f\t%2i\r\n',[left_electrode_x,left_electrode_y,left_electrode_z,left_electrode_species]');
	fprintf(fileID,'\n%s-atom right electrode with lattice param: %s\t%s\t%s\t90.\t90.\t90.\n',...
		num2str(electrode_num_atoms+electrode_dopant_atoms),num2str(cell_size_x),num2str(cell_size_y),num2str(last_right_C_atom_coord-first_right_C_atom_coord+some_space));	
	fprintf(fileID,'%4.5f\t%4.5f\t%4.5f\t%2i\r\n',[right_electrode_x,right_electrode_y,right_electrode_z,right_electrode_species]');
	fprintf(fileID,'\n%s-atom full model (left+device+right) with lattice param: %s\t%s\t%s\t90.\t90.\t90.\n',num2str(num_of_atoms),num2str(cell_size_x),...
		num2str(cell_size_y),num2str(cell_size_z));
	fprintf(fileID,'%4.5f\t%4.5f\t%4.5f\t%2i\r\n',[full_model_x,full_model_y,full_model_z,full_model_species]');
	fclose(fileID);
	
	case 4			% band structure calcs
		% Sorting C atoms first
		C_index = find(species_num == 4);
		x_C = x(C_index);
		y_C = y(C_index);
		z_C = z(C_index);
		len_x = max(x_C) - min(x_C);			% circumferential length
		len_y = max(y_C) - min(y_C);
		len_z = max(z_C) - min(z_C);			% axial length
		
		if chirality(1) == chirality(2)
			some_space_x = d_CC;
			some_space_z = some_space;
		else
			some_space_x = sqrt(3)/2*d_CC;
			some_space_z = some_space;
		end
		
		cell_size_x = len_x + some_space_x;			% 
		cell_size_y = len_y + 25;
		cell_size_z = len_z + some_space;	% 
		
		x = x - x_C(1);
		z = z - z_C(1);
		
		% Sorting dopant atoms
		dopant_index = find(species_num ~= 4);		% find only the dopant molecule atoms

		x_sign = sign(x(dopant_index));
		y_sign = sign(y(dopant_index));
		z_sign = sign(z(dopant_index));
		if sum(x_sign) ~= length(x_sign)*1
			negative_x = find(x_sign<0);
			x(dopant_index(negative_x)) = cell_size_x + x(dopant_index(negative_x));
		end
		
		if sum(y_sign) ~= length(y_sign)*1
			negative_y = find(y_sign<0);
			y(dopant_index(negative_y)) = cell_size_y + y(dopant_index(negative_y));
		end
		
		if sum(z_sign) ~= length(z_sign)*1
			negative_z = find(z_sign<0);
			z(dopant_index(negative_z)) = cell_size_z + z(dopant_index(negative_z));
		end
		
		lat_const = cell_size_z;
		
		lattice_vec_x = [cell_size_x/lat_const 0 0];
		lattice_vec_y = [0 cell_size_y/lat_const 0];
		lattice_vec_z = [0 0 cell_size_z/lat_const];
		
		fileID = fopen(file_to_write,'w');
		fprintf(fileID,'%s\n\n',['from ',file_to_read,'; num of atoms = ',num2str(num_of_atoms),'; lattice const = ',num2str(lat_const)]);
		fprintf(fileID,'%s\n%3.4f\t%3.4f\t%3.4f\n%3.4f\t%3.4f\t%3.4f\n%3.4f\t%3.4f\t%3.4f\n\n','lattice vecs:',[lattice_vec_x;lattice_vec_y;lattice_vec_z]);
		fprintf(fileID,'%s\t\t%s\t\t%s\t\t%s\n','x','y','z','chemical species number');
		fprintf(fileID,'%4.5f\t%4.6f\t%4.6f\t %2i\r\n',[x,y,z,species_num]');
		fclose(fileID);
end
