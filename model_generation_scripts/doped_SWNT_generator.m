% This code generates tube geometry with dopant
% Units used are in angstroms
% The maximum distance of two carbons (maximum width) in the hexagonal 
% lattice is 2*d_CC, the minimum is sqrt(3)/2 * d_CC
% i.e opposite carbon atoms
clear

%% Create CNT
xBox = 45;						% for box generation
yBox = 45;
box = [xBox;yBox];
d_CC  = 1.421;					% bond length in graphite -- obtained from http://phycomp.technion.ac.il/~david/thesis/node3.html
unit = d_CC * 3;
l    = unit * 6;

chirality = [5,5];
CNT = nanotube('C', chirality, d_CC, l);      % generates tube object
pt1 = [+12.769799000  +27.528553000];		% what are these numbers and where do they come from?
pt2 = [+15.000000000  +26.604773000];
CNT = zrotate(CNT,pt1,pt2);        % fixing the points in x direction
CNT = xyMove(CNT, xBox/2, yBox/2); % tube to centre of box

% Rotate CNT tube for armchair configuration so that the hexagons face the dopants tangently
unused_1 = 0;	% unused argument value; some functions require values that are irrelevant in this script
num_of_SWNTs = 1;
CNT = CNT_rotator(num_of_SWNTs,chirality,box,unused_1,CNT);

%% Select dopants
dopant_type = 'AuBr4';
d_dopants = 3;							% number of hexagonal graphite sublattices separating the dopants; 1 = next sublattice, 2 = skip one sublattice...
num_dopants = 4;						% total number of dopants
output = dopant_selector(dopant_type);	% returns a struct containing the reference coords and dopant_gen function number

ref_coords = output.coords;
dopant_gen_fcn = output.fcn_num;

%% Create dopant_1
offset_1 = 1;					% the offset count starts from z = 0, and goes by number of unit cell, i.e. offset = 2 skips the first unit cell; check the center_coord loop for details
position_1 = 1;					% 1 = top of CNT; -1 = bottom of CNT

% Set orientation of dopant molecule
x_rot_1 = -90;
y_rot_1 = 45+90;
z_rot_1 = 0;
dopant_dist_1 = 5;				% center to center dist in angstroms

% Generating/duplicating the 'center' atom (reference atom) coordinates of the dopant
center_coords_1 = ref_atom_coords_gen(chirality,box,CNT,position_1,num_dopants,d_dopants,offset_1,dopant_dist_1);

% Rotating the dopant with the 'center' atom (reference atom) as the origin
rot_matrix_1 = dopant_orientation(x_rot_1,y_rot_1,z_rot_1);
new_ref_coords_1 = (rot_matrix_1*ref_coords')';

% Duplicating the rest of the atoms in the dopant molecule
if dopant_gen_fcn == 1
	dopant_1 = KAuBr4_gen(new_ref_coords_1,center_coords_1);
elseif dopant_gen_fcn == 2
	dopant_1 = AuBr4_gen(new_ref_coords_1,center_coords_1);
else
	dopant_1 = K_gen(new_ref_coords_1,center_coords_1);
end

atoms_arg_1 = num2cell(dopant_1.atoms);

%% Create dopant_2
offset_2 = 1;					% the offset count starts from z = 0, and goes by number of unit cell, i.e. offset = 2 skips the first unit cell; check the center_coord loop for details
position_2 = -1;					% 1 = top of CNT; -1 = bottom of CNT

% Set orientation of dopant molecule
x_rot_2 = -90;
y_rot_2 = 45+90;
z_rot_2 = 180;
dopant_dist_2 = 5;				% center to center dist in angstroms

% Generating/duplicating the 'center' atom (reference atom) coordinates of the dopant
center_coords_2 = ref_atom_coords_gen(chirality,box,CNT,position_2,num_dopants,d_dopants,offset_2,dopant_dist_2);

% Rotating the dopant with the 'center' atom (reference atom) as the origin
rot_matrix_2 = dopant_orientation(x_rot_2,y_rot_2,z_rot_2);
new_ref_coords_2 = (rot_matrix_2*ref_coords')';

% Duplicating the rest of the atoms in the dopant molecule
if dopant_gen_fcn == 1
	dopant_2 = KAuBr4_gen(new_ref_coords_2,center_coords_2);
elseif dopant_gen_fcn == 2
	dopant_2 = AuBr4_gen(new_ref_coords_2,center_coords_2);
else
	dopant_2 = K_gen(new_ref_coords_2,center_coords_2);
end

atoms_arg_2 = num2cell(dopant_2.atoms);

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

outxyz([folder_name,'/',filename],CNT,atoms_arg_1{:},atoms_arg_2{:}); 
cmd = sprintf('xcrysden --xyz coord/%s/%s.xyz &',folder_name,filename);
system(cmd);