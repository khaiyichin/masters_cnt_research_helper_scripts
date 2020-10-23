% This code generates junction (or 2 tubes) geometry with dopant
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
l    = unit * 9;

intertube = 6.6883;				% distance between the two tubes
overlap = l;					% overlap of tubes at junction; use unit * num for incomplete overlap

chirality = [5,5];
CNT = nanotube('C', chirality, d_CC, l);      % generates tube object
pt1 = [+12.769799000  +27.528553000];
pt2 = [+15.000000000  +26.604773000];
CNT = zrotate(CNT,pt1,pt2);        %fixing the points in x direction
CNT = xyMove(CNT, xBox/2, yBox/2); %tube to centre of box

%% Create junction from CNT
junc = junction(CNT.coord);	% junction is a user-defined class
junc.overLap = overlap;
junc.unit = unit;
junc.intTube = intertube;
junc.sysL = l;
junc.symbol = 'C';
junc.diameter = CNT.diameter;
junc.upperOffset = 0; %1.421 * 3^0.5 / 2;
junc = MakeJunc(junc);

% Rotate CNT tube for armchair configuration so that the hexagons face the dopants tangently
unused_1 = 0;	% unused argument value; some functions require values that are irrelevant in this script
num_of_SWNTs = 2;
junc = CNT_rotator(num_of_SWNTs,chirality,box,intertube,junc);

%% Select dopants
dopant_type = 'KAuBr4';
d_dopants = 3;							% number of hexagonal graphite sublattices separating the dopants; 1 = next sublattice, 2 = skip one sublattice...
num_dopants = 4;						% total number of dopants
output = dopant_selector(dopant_type);	% returns a struct containing the reference coords and dopant_gen function number

ref_coords = output.coords;
dopant_gen_fcn = output.fcn_num;

% ref_coords = ...				% temp (8,0) KAuBr4 dopant
% 	[	0         0         0
% 		0.0008   -0.5414   -4.2703
% 		-1.7837    0.0072   -1.7460
% 		-1.7922   -0.0199    1.7411
% 		1.7926   -0.0178    1.7411
% 		1.7834    0.0100   -1.7464];

% ref_coords = ...				% temp (5,5) KAuBr4 dopant
%     [	0         0         0
% 		-0.0005    1.3740   -3.6971
% 		-1.8379   -0.0050   -1.7276
% 		-1.8439    0.0010    1.7218
% 	    1.8440    0.0026    1.7221
% 		1.8381   -0.0027   -1.7274];

%% Create dopant_1
offset_1 = 2;						% the offset count starts from z = 0, and goes by number of unit cell, i.e. offset = 2 skips the first unit cell; check the center_coord loop for details

% Set orientation of dopant molecule
x_rot_1 = 90;
y_rot_1 = 45;
z_rot_1 = 180;

% Set translation parameters (in angstroms)
x_translation_1 = 0;
y_translation_1 = 0;
z_translation_1 = 0;

% Generating/duplicating the 'center' atom (reference atom) coordinates of the dopant
center_coords_1 = ref_atom_coords_gen(chirality,box,CNT,unused_1,num_dopants,d_dopants,offset_1,unused_1);
array_size = size(center_coords_1);
translation_1 = dopant_translation(x_translation_1,y_translation_1,z_translation_1,array_size);

% Rotating the dopant with the 'center' atom (reference atom) as the origin
rot_matrix_1 = dopant_orientation(x_rot_1,y_rot_1,z_rot_1);
new_ref_coords_1 = (rot_matrix_1*ref_coords')';

% Set translation of dopant molecule duplication
center_coords_1 = center_coords_1 + translation_1;

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
offset_2 = 2;						% the offset count starts from z = 0, and goes by number of unit cell, i.e. offset = 2 skips the first unit cell; check the center_coord loop for details

% Set orientation of dopant molecule
x_rot_2 = 90;
y_rot_2 = 45;
z_rot_2 = 180;

% Set translation parameters (in angstroms)
x_translation_2 = 0;
y_translation_2 = 13.47305787;
z_translation_2 = 0;

% Generating/duplicating the 'center' atom (reference atom) coordinates of the dopant
center_coords_2 = ref_atom_coords_gen(chirality,box,CNT,unused_1,num_dopants,d_dopants,offset_2,unused_1);
array_size = size(center_coords_2);
translation_2 = dopant_translation(x_translation_2,y_translation_2,z_translation_2,array_size);

% Rotating the dopant with the 'center' atom (reference atom) as the origin
rot_matrix_2 = dopant_orientation(x_rot_2,y_rot_2,z_rot_2);
new_ref_coords_2 = (rot_matrix_2*ref_coords')';

% Set translation of dopant molecule duplication
center_coords_2 = center_coords_2 + translation_2;

% Duplicating the rest of the atoms in the dopant molecule
if dopant_gen_fcn == 1
	dopant_2 = KAuBr4_gen(new_ref_coords_2,center_coords_2);
elseif dopant_gen_fcn == 2
	dopant_2 = AuBr4_gen(new_ref_coords_2,center_coords_2);
else
	dopant_2 = K_gen(new_ref_coords_2,center_coords_2);
end

atoms_arg_2 = num2cell(dopant_2.atoms);

%% Create dopant_3
offset_3 = 2;						% the offset count starts from z = 0, and goes by number of unit cell, i.e. offset = 2 skips the first unit cell; check the center_coord loop for details

% Set orientation of dopant molecule
x_rot_3 = 90;
y_rot_3 = 45;
z_rot_3 = 0;

% Set translation parameters (in angstroms)
x_translation_3 = 0;
y_translation_3 = -13.47305787;
z_translation_3 = 0;

% Generating/duplicating the 'center' atom (reference atom) coordinates of the dopant
center_coords_3 = ref_atom_coords_gen(chirality,box,CNT,unused_1,num_dopants,d_dopants,offset_3,unused_1);
array_size = size(center_coords_3);
translation_3 = dopant_translation(x_translation_3,y_translation_3,z_translation_3,array_size);

% Rotating the dopant with the 'center' atom (reference atom) as the origin
rot_matrix_3 = dopant_orientation(x_rot_3,y_rot_3,z_rot_3);
new_ref_coords_3 = (rot_matrix_3*ref_coords')';

% Set translation of dopant molecule duplication
center_coords_3 = center_coords_3 + translation_3;

% Duplicating the rest of the atoms in the dopant molecule
if dopant_gen_fcn == 1
	dopant_3 = KAuBr4_gen(new_ref_coords_3,center_coords_3);
elseif dopant_gen_fcn == 2
	dopant_3 = AuBr4_gen(new_ref_coords_3,center_coords_3);
else
	dopant_3 = K_gen(new_ref_coords_3,center_coords_3);
end

atoms_arg_3 = num2cell(dopant_3.atoms);

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
filename = [hour,':',min,':',sec,'_',chirality,'junc'];
if ~exist(['coord/',folder_name], 'dir')
  mkdir(['coord/',folder_name]);
end

outxyz([folder_name,'/',filename],junc,atoms_arg_1{:},atoms_arg_2{:},atoms_arg_3{:})
cmd = sprintf('xcrysden --xyz coord/%s/%s.xyz &',folder_name,filename);
system(cmd);
