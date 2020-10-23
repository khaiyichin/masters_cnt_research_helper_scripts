function CNT_object = CNT_rotator(num_of_tubes,chirality,box,intertube,CNT_object)
%CNT_rotator Rotates SWNT(s) for armchair configuration so that the 
%hexagons face the dopants tangently
%   This works for the (5,5) armchair configuration (theoretically all 
%	armchair config), but for anything else you must check to see if it's 
%	still right.

d_CC = 1.421;
xBox = box(1);
yBox = box(2);

if num_of_tubes == 1
		if chirality(1) == chirality(2)
		num_atoms_per_layer_per_tube = chirality(1) + chirality(2);

		tube_1st_layer = CNT_object.coord(1:num_atoms_per_layer_per_tube,:);

% 		diam = 3*chirality(1)/pi*d_CC;			% diameter of tube -- obtained from http://www.photon.t.u-tokyo.ac.jp/~maruyama/kataura/chirality.html
		center = [xBox/2,yBox/2,0];		% center of tube

		tube_trans = tube_1st_layer - center;
		angles = atan2(tube_trans(:,2),tube_trans(:,1));

		for i = 1:length(angles)
			if angles(i) < 0
				angles(i) = angles(i) + 2*pi;
			end
		end

		% Taking angles from the 90 degree line (+ y-axis)
		[~,ind_1] = min(abs(angles-pi/2));
		ang_1 = angles(ind_1);
		angles(ind_1) = NaN;		% remove this angle for consideration of minimum in the next few lines

		[~,ind_2] = min(abs(angles-pi/2));
		ang_2 = angles(ind_2);

		while norm(tube_1st_layer(ind_1,:)-tube_1st_layer(ind_2,:),2) > d_CC	% check to see if the two points (atoms) found are directly connected to each other
			angles(ind_2) = NaN;
			[~,ind_2] = min(abs(angles-pi/2));
			ang_2 = angles(ind_2);
		end

		ang_between_pt_1_2 = abs(ang_2 - ang_1);
		dir_to_rotate = sign(pi/2-ang_1);
		ang_to_rotate = dir_to_rotate*(abs(pi/2-ang_1) + ang_between_pt_1_2/2);

		% Updating coordinates for CNT
		transformed_tube = CNT_object.coord-center;
		rot_z_correction = [cos(ang_to_rotate) -sin(ang_to_rotate) 0; sin(ang_to_rotate) cos(ang_to_rotate) 0; 0 0 1];
		new_coord = (rot_z_correction*transformed_tube')';
		new_coord = new_coord(1:length(transformed_tube(:,1)),:) + center;
		CNT_object.coord = coord_zsort(new_coord);
	end
	
elseif num_of_tubes == 2
		if chirality(1) == chirality(2)
		num_atoms_per_layer_per_tube = chirality(1) + chirality(2);

		first_tube_1st_layer = CNT_object.coord(1:num_atoms_per_layer_per_tube,:);
		second_tube_1st_layer = CNT_object.coord(num_atoms_per_layer_per_tube+1:2*num_atoms_per_layer_per_tube,:);

		diam = 3*chirality(1)/pi*d_CC;			% diameter of tube -- obtained from http://www.photon.t.u-tokyo.ac.jp/~maruyama/kataura/chirality.html
		center_first = [xBox/2,yBox/2-intertube/2-diam/2,0];		% center of first tube
		center_second = [xBox/2,yBox/2+intertube/2+diam/2,0];		% center of second tube

		% FIRST TUBE (BOTTOM)
		first_tube_trans = first_tube_1st_layer - center_first;
		angles_first = atan2(first_tube_trans(:,2),first_tube_trans(:,1));

		for i = 1:length(angles_first)
			if angles_first(i) < 0
				angles_first(i) = angles_first(i) + 2*pi;
			end
		end

		% Taking angles from the 90 degree line (+ y-axis)
		[~,ind_1_first] = min(abs(angles_first-pi/2));
		ang_1_first = angles_first(ind_1_first);
		angles_first(ind_1_first) = NaN;		% remove this angle for consideration of minimum in the next few lines

		[~,ind_2_first] = min(abs(angles_first-pi/2));
		ang_2_first = angles_first(ind_2_first);

		while norm(first_tube_1st_layer(ind_1_first,:)-first_tube_1st_layer(ind_2_first,:),2) > d_CC	% check to see if the two points (atoms) found are directly connected to each other
			angles_first(ind_2_first) = NaN;
			[~,ind_2_first] = min(abs(angles_first-pi/2));
			ang_2_first = angles_first(ind_2_first);
		end

		ang_between_pt_1_2_first = abs(ang_2_first - ang_1_first);
		dir_to_rotate_first = sign(pi/2-ang_1_first);
		ang_to_rotate_first = dir_to_rotate_first*(abs(pi/2-ang_1_first) + ang_between_pt_1_2_first/2);

		% SECOND TUBE (TOP)
		second_tube_trans = second_tube_1st_layer - center_second;
		angles_second = atan2(second_tube_trans(:,2),second_tube_trans(:,1));

		for i = 1:length(angles_second)
			if angles_second(i) < 0
				angles_second(i) = angles_second(i) + 2*pi;
			end
		end

		% Taking angles from the 270 degree line (- y-axis)
		[~,ind_1_second] = min(abs(angles_second-3*pi/2));
		ang_1_second = angles_second(ind_1_second);
		angles_second(ind_1_second) = NaN;		% remove this angle for consideration of minimum in the next few lines

		[~,ind_2_second] = min(abs(angles_second-3*pi/2));
		ang_2_second = angles_second(ind_2_second);

		while norm(second_tube_1st_layer(ind_1_second,:)-second_tube_1st_layer(ind_2_second,:),2) > d_CC	% check to see if the two points (atoms) found are directly connected to each other
			angles_second(ind_2_second) = NaN;
			[~,ind_2_second] = min(abs(angles_second-3*pi/2));
			ang_2_second = angles_first(ind_2_second);
		end

		ang_between_pt_1_2_second = abs(ang_2_second - ang_1_second);
		dir_to_rotate_second = sign(3*pi/2-ang_1_second);
		ang_to_rotate_second =  dir_to_rotate_second*(abs(3*pi/2-ang_1_second) + ang_between_pt_1_2_second/2);

		% Updating coordinates for both CNTs
		first_tube = first_tube_1st_layer;
		second_tube = second_tube_1st_layer;

		layers = length(CNT_object.coord(:,1))/num_atoms_per_layer_per_tube;

		for i = 1:layers/2 -1
			first_tube = [first_tube;CNT_object.coord((2*i*num_atoms_per_layer_per_tube+1):(2*i+1)*num_atoms_per_layer_per_tube,:)];
			second_tube = [second_tube;CNT_object.coord((2*i+1)*num_atoms_per_layer_per_tube+1:(2*i+2)*num_atoms_per_layer_per_tube,:)];
		end

		new_first_tube = first_tube - center_first;
		new_second_tube = second_tube - center_second;

		rot_z_first = [cos(ang_to_rotate_first) -sin(ang_to_rotate_first) 0; sin(ang_to_rotate_first) cos(ang_to_rotate_first) 0; 0 0 1];
		rot_z_second = [cos(ang_to_rotate_second) -sin(ang_to_rotate_second) 0; sin(ang_to_rotate_second) cos(ang_to_rotate_second) 0; 0 0 1];
		new_coord = [rot_z_first*new_first_tube',rot_z_second*new_second_tube']';
		new_coord = [new_coord(1:length(new_first_tube(:,1)),:) + center_first ; new_coord(length(new_first_tube(:,1))+1:end,:) + center_second];
		CNT_object.coord = coord_zsort(new_coord);

	end
	
end

end

