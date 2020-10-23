% This code generates tube geometry with dopant
% Units used are in angstroms
% The maximum distance of two carbons (maximum width) in the hexagonal 
% lattice is 2*d_CC, the minimum is sqrt(3)/2 * d_CC
% i.e opposite carbon atoms
clear

%% Create CNT
xBox = 45;						% for box generation
yBox = 45;
box_center = [xBox/2, yBox/2, 0];
d_CC  = 1.421;					% bond length in graphite -- obtained from http://phycomp.technion.ac.il/~david/thesis/node3.html
unit = d_CC * 3;
l    = unit * 6;

d_dopants = 1;					% number of hexagonal graphite sublattices separating the dopants; 1 = next sublattice, 2 = skip one sublattice...
num_dopants = 10;				% total number of dopants
siestaOutDir = ['./siesta']; %in local siesta directory

chirality = [5,5];
CNT = nanotube('C', chirality, d_CC, l);      % generates tube object
pt1 = [+12.769799000  +27.528553000];		% what are these numbers and where do they come from?
pt2 = [+15.000000000  +26.604773000];
CNT = zrotate(CNT,pt1,pt2);        %fixing the points in x direction
CNT = xyMove(CNT, xBox/2, yBox/2); %tube to centre of box

%% Rotate CNT tube for armchair configuration so that the hexagons face the dopants tangently
% This works for the (5,5) armchair configuration (theoretically all armchair config)
% , but for anything else you must check to see if it's still right
if chirality(1) == chirality(2)
	num_atoms_per_layer_per_tube = chirality(1) + chirality(2);
	
	tube_1st_layer = CNT.coord(1:num_atoms_per_layer_per_tube,:);
	
	diam = 3*chirality(1)/pi*d_CC;			% diameter of tube -- obtained from http://www.photon.t.u-tokyo.ac.jp/~maruyama/kataura/chirality.html
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
	transformed_tube = CNT.coord-center;
	rot_z_correction = [cos(ang_to_rotate) -sin(ang_to_rotate) 0; sin(ang_to_rotate) cos(ang_to_rotate) 0; 0 0 1];
	new_coord = (rot_z_correction*transformed_tube')';
	new_coord = new_coord(1:length(transformed_tube(:,1)),:) + center;
	CNT.coord = coord_zsort(new_coord);
end

%% Create KAuBr4 dopant
% ref_coords = ...				% KAuBr4 dopant
% 	[	0.0000	0.0000	0.0000	;
% 		2.1616	-2.1373	2.7727	;
% 		-0.0134	-2.4093	0.0000	;
% 		2.4146	0.0000	0.0000	;
% 		0.0134	2.4093	0.0000	;
% 		-2.4146	0.0000	0.0000	];
	
% ref_coords = ...				% AuBr4 dopant
% 	[	0.0000	0.0000	0.0000	;
% 		-0.0134	-2.4093	0.0000	;
% 		2.4146	0.0000	0.0000	;
% 		0.0134	2.4093	0.0000	;
% 		-2.4146	0.0000	0.0000	]; 

ref_coords = ...				% K dopant
	[	0.000	0.000	0.000	];

offset = 1;						% the offset count starts from z = 0, and goes by number of unit cell, i.e. offset = 2 skips the first unit cell; check the center_coord loop for details
potassium = dopant;				% dopant is a user-defined class
potassium.symbol = 'K';
d_K_to_CNT = 1;					% center to center dist, [Angstroms]

gold = dopant;
gold.symbol = 'Au';
d_gold_to_CNT = 1;				% center to center dist, [Angstroms]

bromine_1 = dopant;
bromine_1.symbol = 'Br';

bromine_2 = dopant;
bromine_2.symbol = 'Br';

bromine_3 = dopant;
bromine_3.symbol = 'Br';

bromine_4 = dopant;
bromine_4.symbol = 'Br';

% Now you only have to duplicate the KAuBr4 (the center atom of KAuBr4) and
% the rest will duplicate itself
% KAuBr4 = molecule([gold,potassium,bromine_1,bromine_2,bromine_3,bromine_4],gold);
% AuBr4 = molecule([gold,bromine_1,bromine_2,bromine_3,bromine_4],gold);
K = molecule(potassium,potassium);

% Generate coordinates for center atom of dopant -- strictly for armchair
% or zig-zag configuration only; for other configuration these coordinates
% might be incorrect
if chirality(1) == chirality(2) % armchair config
	x_coord = xBox/2;
	y_coord = yBox/2 + CNT.diameter/2 + d_gold_to_CNT;
	z_coord = sqrt(3)/2*d_CC;
	
	center_coords = [];
	for i = offset:d_dopants:(num_dopants*d_dopants)-(d_dopants-offset)
		center_coords = [	center_coords;
							x_coord, y_coord, (2*i-1)*z_coord	];
	end
	
else % zig-zag config
	rot_z_correction = eye(3);
	x_coord = xBox/2;
	y_coord = yBox/2 + CNT.diameter/2 + d_gold_to_CNT;
	z_coord = d_CC;
	
	center_coords = [];
	for i = offset:d_dopants:(num_dopants*d_dopants)-(d_dopants-offset)
		center_coords = [	center_coords;
							x_coord, y_coord, (3*i-2)*z_coord	];
	end
	
end

% Set orientation of dopant molecule
rot_ang_x = deg2rad(-90);		% the CNT axis is along z, whereas the dopant is square planar in x-y
rot_x = [1 0 0; 0 cos(rot_ang_x) -sin(rot_ang_x); 0 sin(rot_ang_x) cos(rot_ang_x)];
rot_ang_y = deg2rad(45+90);
rot_y = [cos(rot_ang_y) 0 sin(rot_ang_y); 0 1 0; -sin(rot_ang_y) 0 cos(rot_ang_y)];
rot_ang_z = deg2rad(0);
rot_z = [cos(rot_ang_z) -sin(rot_ang_z) 0; sin(rot_ang_z) cos(rot_ang_z) 0; 0 0 1];
new_ref_coords = (rot_z*rot_y*rot_x*ref_coords')';

% KAuBr4.relative_pos(new_ref_coords);
% KAuBr4.center_atom_duplication(center_coords);
% KAuBr4.rem_atoms_duplication();
% AuBr4.relative_pos(new_ref_coords);
% AuBr4.center_atom_duplication(center_coords);
% AuBr4.rem_atoms_duplication();
K.relative_pos(ref_coords);
K.center_atom_duplication(center_coords);

% atoms_arg = num2cell(KAuBr4.atoms);
% atoms_arg = num2cell(AuBr4.atoms);
atoms_arg = num2cell(K.atoms);

%% Create alternate_1 KAuBr4 molecules
offset_alt_1 = 2;						% the offset count starts from z = 0, and goes by number of unit cell, i.e. offset = 2 skips the first unit cell; check the center_coord loop for details

potassium_alt_1 = dopant;				% dopant is a user-defined class
potassium_alt_1.symbol = 'K';

gold_alt_1 = dopant;
gold_alt_1.symbol = 'Au';

bromine_1_alt_1 = dopant;
bromine_1_alt_1.symbol = 'Br';

bromine_2_alt_1 = dopant;
bromine_2_alt_1.symbol = 'Br';

bromine_3_alt_1 = dopant;
bromine_3_alt_1.symbol = 'Br';

bromine_4_alt_1 = dopant;
bromine_4_alt_1.symbol = 'Br';

% Now you only have to duplicate the KAuBr4 (the center atom of KAuBr4) and
% the rest will duplicate itself
% KAuBr4_alt_1 = molecule([gold_alt_1,potassium_alt_1,bromine_1_alt_1,bromine_2_alt_1,bromine_3_alt_1,bromine_4_alt_1],gold_alt_1);
% AuBr4_alt_1 = molecule([gold_alt_1,bromine_1_alt_1,bromine_2_alt_1,bromine_3_alt_1,bromine_4_alt_1],gold_alt_1);


% Generate coordinates for center atom of dopant -- strictly for armchair
% or zig-zag configuration only; for other configuration these coordinates
% might be incorrect
if chirality(1) == chirality(2) % armchair config
	x_coord = xBox/2;
	y_coord = yBox/2 + CNT.diameter/2 + d_gold_to_CNT;
	z_coord = sqrt(3)/2*d_CC;
	
	center_coords_alt_1 = [];
	for i = offset_alt_1:d_dopants:(num_dopants*d_dopants)-(d_dopants-offset_alt_1)
		center_coords_alt_1 = [	center_coords_alt_1;
							x_coord, y_coord, (2*i-1)*z_coord	];
	end
	
else % zig-zag config
	x_coord = xBox/2;
	y_coord = yBox/2 + CNT.diameter/2 + d_gold_to_CNT;
	z_coord = d_CC;
	
	center_coords_alt_1 = [];
	for i = offset_alt_1:d_dopants:(num_dopants*d_dopants)-(d_dopants-offset_alt_1)
		center_coords_alt_1 = [	center_coords_alt_1;
							x_coord, y_coord, (3*i-2)*z_coord	];
	end
	
end

% Set orientation of dopant molecule
rot_ang_x_alt_1 = deg2rad(-90);		% the CNT axis is along z, whereas the dopant is square planar in x-y
rot_x_alt_1 = [1 0 0; 0 cos(rot_ang_x_alt_1) -sin(rot_ang_x_alt_1); 0 sin(rot_ang_x_alt_1) cos(rot_ang_x_alt_1)];
rot_ang_y_alt_1 = deg2rad(45+270);
rot_y_alt_1 = [cos(rot_ang_y_alt_1) 0 sin(rot_ang_y_alt_1); 0 1 0; -sin(rot_ang_y_alt_1) 0 cos(rot_ang_y_alt_1)];
rot_ang_z_alt_1 = deg2rad(0);
rot_z_alt_1 = [cos(rot_ang_z_alt_1) -sin(rot_ang_z_alt_1) 0; sin(rot_ang_z_alt_1) cos(rot_ang_z_alt_1) 0; 0 0 1];
new_ref_coords_alt_1 = (rot_z_alt_1*rot_y_alt_1*rot_x_alt_1*ref_coords')';

% KAuBr4_alt_1.relative_pos(new_ref_coords_alt_1);
% KAuBr4_alt_1.center_atom_duplication(center_coords_alt_1);
% KAuBr4_alt_1.rem_atoms_duplication();
% AuBr4_alt_1.relative_pos(new_ref_coords_alt_1);
% AuBr4_alt_1.center_atom_duplication(center_coords_alt_1);
% AuBr4_alt_1.rem_atoms_duplication();

% atoms_arg_alt_1 = num2cell(KAuBr4_alt_1.atoms);
% atoms_arg_alt_1 = num2cell(AuBr4_alt_1.atoms);

%% Create alternate_2 KAuBr4 molecules
offset_alt_2 = 2;						% the offset count starts from z = 0, and goes by number of unit cell, i.e. offset = 2 skips the first unit cell; check the center_coord loop for details

potassium_alt_2 = dopant;				% dopant is a user-defined class
potassium_alt_2.symbol = 'K';

gold_alt_2 = dopant;
gold_alt_2.symbol = 'Au';

bromine_1_alt_2 = dopant;
bromine_1_alt_2.symbol = 'Br';

bromine_2_alt_2 = dopant;
bromine_2_alt_2.symbol = 'Br';

bromine_3_alt_2 = dopant;
bromine_3_alt_2.symbol = 'Br';

bromine_4_alt_2 = dopant;
bromine_4_alt_2.symbol = 'Br';

% Now you only have to duplicate the KAuBr4 (the center atom of KAuBr4) and
% the rest will duplicate itself
% KAuBr4_alt_2 = molecule([gold_alt_2,potassium_alt_2,bromine_1_alt_2,bromine_2_alt_2,bromine_3_alt_2,bromine_4_alt_2],gold_alt_2);
% AuBr4_alt_2 = molecule([gold_alt_2,bromine_1_alt_2,bromine_2_alt_2,bromine_3_alt_2,bromine_4_alt_2],gold_alt_2);

% Generate coordinates for center atom of dopant -- strictly for armchair
% or zig-zag configuration only; for other configuration these coordinates
% might be incorrect
if chirality(1) == chirality(2) % armchair config
	x_coord = xBox/2;
	y_coord = yBox/2 - CNT.diameter/2 - d_gold_to_CNT;
	z_coord = sqrt(3)*d_CC;
	
	center_coords_alt_2 = [];
	for i = offset_alt_2:d_dopants:(num_dopants*d_dopants)-(d_dopants-offset_alt_2)
		center_coords_alt_2 = [	center_coords_alt_2;
							x_coord, y_coord, i*z_coord	];
	end
	
else % zig-zag config
	x_coord = xBox/2;
	y_coord = yBox/2 - CNT.diameter/2 - d_gold_to_CNT;
	z_coord = d_CC;
	
	center_coords_alt_2 = [];
	for i = offset_alt_2:d_dopants:(num_dopants*d_dopants)-(d_dopants-offset_alt_2)
		center_coords_alt_2 = [	center_coords_alt_2;
							x_coord, y_coord, (3*i-2)*z_coord	];
	end
	
end

% Set orientation of dopant molecule
rot_ang_x_alt_2 = deg2rad(-90);		% the CNT axis is along z, whereas the dopant is square planar in x-y
rot_x_alt_2 = [1 0 0; 0 cos(rot_ang_x_alt_2) -sin(rot_ang_x_alt_2); 0 sin(rot_ang_x_alt_2) cos(rot_ang_x_alt_2)];
rot_ang_y_alt_2 = deg2rad(45+90);
rot_y_alt_2 = [cos(rot_ang_y_alt_2) 0 sin(rot_ang_y_alt_2); 0 1 0; -sin(rot_ang_y_alt_2) 0 cos(rot_ang_y_alt_2)];
rot_ang_z_alt_2 = deg2rad(180);
rot_z_alt_2 = [cos(rot_ang_z_alt_2) -sin(rot_ang_z_alt_2) 0; sin(rot_ang_z_alt_2) cos(rot_ang_z_alt_2) 0; 0 0 1];
new_ref_coords_alt_2 = (rot_z_alt_2*rot_y_alt_2*rot_x_alt_2*ref_coords')';

% KAuBr4_alt_2.relative_pos(new_ref_coords_alt_2);
% KAuBr4_alt_2.center_atom_duplication(center_coords_alt_2);
% KAuBr4_alt_2.rem_atoms_duplication();
% AuBr4_alt_2.relative_pos(new_ref_coords_alt_2);
% AuBr4_alt_2.center_atom_duplication(center_coords_alt_2);
% AuBr4_alt_2.rem_atoms_duplication();

% atoms_arg_alt_2 = num2cell(KAuBr4_alt_2.atoms);
% atoms_arg_alt_2 = num2cell(AuBr4_alt_2.atoms);

%% Create alternate_2 KAuBr4 molecules
offset_alt_3 = 7;						% the offset count starts from z = 0, and goes by number of unit cell, i.e. offset = 2 skips the first unit cell; check the center_coord loop for details

potassium_alt_3 = dopant;				% dopant is a user-defined class
potassium_alt_3.symbol = 'K';

gold_alt_3 = dopant;
gold_alt_3.symbol = 'Au';

bromine_1_alt_3 = dopant;
bromine_1_alt_3.symbol = 'Br';

bromine_2_alt_3 = dopant;
bromine_2_alt_3.symbol = 'Br';

bromine_3_alt_3 = dopant;
bromine_3_alt_3.symbol = 'Br';

bromine_4_alt_3 = dopant;
bromine_4_alt_3.symbol = 'Br';

% Now you only have to duplicate the KAuBr4 (the center atom of KAuBr4) and
% the rest will duplicate itself
% KAuBr4_alt_3 = molecule([gold_alt_3,potassium_alt_3,bromine_1_alt_3,bromine_2_alt_3,bromine_3_alt_3,bromine_4_alt_3],gold_alt_3);
% AuBr4_alt_3 = molecule([gold_alt_3,bromine_1_alt_3,bromine_2_alt_3,bromine_3_alt_3,bromine_4_alt_3],gold_alt_3);

% Generate coordinates for center atom of dopant -- strictly for armchair
% or zig-zag configuration only; for other configuration these coordinates
% might be incorrect
if chirality(1) == chirality(2) % armchair config
	x_coord = xBox/2;
	y_coord = yBox/2 - CNT.diameter/2 - d_gold_to_CNT;
	z_coord = sqrt(3)*d_CC;
	
	center_coords_alt_3 = [];
	for i = offset_alt_3:d_dopants:(num_dopants*d_dopants)-(d_dopants-offset_alt_3)
		center_coords_alt_3 = [	center_coords_alt_3;
							x_coord, y_coord, i*z_coord	];
	end
	
else % zig-zag config
	x_coord = xBox/2;
	y_coord = yBox/2 - CNT.diameter/2 - d_gold_to_CNT;
	z_coord = d_CC;
	
	center_coords_alt_3 = [];
	for i = offset_alt_3:d_dopants:(num_dopants*d_dopants)-(d_dopants-offset_alt_3)
		center_coords_alt_3 = [	center_coords_alt_3;
							x_coord, y_coord, (3*i-2)*z_coord	];
	end
	
end

% Set orientation of dopant molecule
rot_ang_x_alt_3 = deg2rad(-90);		% the CNT axis is along z, whereas the dopant is square planar in x-y
rot_x_alt_3 = [1 0 0; 0 cos(rot_ang_x_alt_3) -sin(rot_ang_x_alt_3); 0 sin(rot_ang_x_alt_3) cos(rot_ang_x_alt_3)];
rot_ang_y_alt_3 = deg2rad(45+270);
rot_y_alt_3 = [cos(rot_ang_y_alt_3) 0 sin(rot_ang_y_alt_3); 0 1 0; -sin(rot_ang_y_alt_3) 0 cos(rot_ang_y_alt_3)];
rot_ang_z_alt_3 = deg2rad(180);
rot_z_alt_3 = [cos(rot_ang_z_alt_3) -sin(rot_ang_z_alt_3) 0; sin(rot_ang_z_alt_3) cos(rot_ang_z_alt_3) 0; 0 0 1];
new_ref_coords_alt_3 = (rot_z_alt_3*rot_y_alt_3*rot_x_alt_3*ref_coords')';

% KAuBr4_alt_3.relative_pos(new_ref_coords_alt_3);
% KAuBr4_alt_3.center_atom_duplication(center_coords_alt_3);
% KAuBr4_alt_3.rem_atoms_duplication();
% AuBr4_alt_3.relative_pos(new_ref_coords_alt_3);
% AuBr4_alt_3.center_atom_duplication(center_coords_alt_3);
% AuBr4_alt_3.rem_atoms_duplication();

% atoms_arg_alt_3 = num2cell(KAuBr4_alt_3.atoms);
% atoms_arg_alt_3 = num2cell(AuBr4_alt_3.atoms);

%% Save to file and open in xcrysden
[hour,min,sec] = hms(datetime('now'));
[year,mon,day] = ymd(datetime('now'));

hour = num2str(hour);
min = num2str(min);
sec = num2str(round(sec));
year = num2str(year);
mon = num2str(mon);
day = num2str(day);
chirality = [num2str(chirality(1)),num2str(chirality(2))];

cd /home/khaiyi/Documents/Data/CNT_with_KAuBr4/

folder_name = [year,mon,day];
filename = [hour,':',min,':',sec,'_',chirality,'tube'];
if ~exist(['coord/',folder_name], 'dir')
  mkdir(['coord/',folder_name]);
end

outxyz([folder_name,'/',filename],CNT,atoms_arg{:});%,atoms_arg_alt_2{:});%,atoms_arg_alt_2{:});%,atoms_arg_alt_3{:});
cmd = sprintf('xcrysden --xyz coord/%s/%s.xyz &',folder_name,filename);
system(cmd);