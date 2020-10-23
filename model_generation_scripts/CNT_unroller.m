% This code 'unwraps' CNTs by generating a graphene sheet.
% Units used are in angstroms
clear

%% Create unrolled CNT
chirality = [5,5];
unit_cell_tube = 3;		% number of unit cells along tube axis (z-axis)
unit_cell_circumference = chirality(1);

d_CC = 1.421;	% bond length in graphite -- obtained from http://phycomp.technion.ac.il/~david/thesis/node3.html

% One graphene unit
if chirality(1) == chirality(2)	% armchair
	unit_x = 3*d_CC;		% zigzag unit cell length
	unit_z = sqrt(3)*d_CC;	% armchair unit cell length
	type = 1;
else	% zigzag
	unit_x = sqrt(3)*d_CC;	% armchair unit cell length
	unit_z = 3*d_CC;		% zigzag unit cell length
	type = 0;
end

% lattice size
xBox = unit_cell_circumference*unit_x; 
zBox = unit_cell_tube*unit_z;
yBox = 0;

% graphene position
origin = [xBox/2,yBox/2,zBox/2];
sheet= unrolled_CNT(unit_cell_circumference,unit_cell_tube,origin,type);

%% Adding dopant from relaxed model
cd /home/khaiyi/Documents/Data/CNT_with_KAuBr4/SWNT/5-5/55SWNT_8KAuBr4_06duc/TACC_calcs/conduction_calcs/lonestar5_calcs/
file_to_read = '55SWNT_8KAuBr4_06duc.xyz';

cmd = sprintf('xcrysden --xyz %s&',file_to_read);
system(cmd);

% Parameters
dopant_type = 'KAuBr4';
x_steps_1 = 1;				% row 1 dopants
z_steps_1 = 2;
x_steps_2 = 6;				% row 2 dopants
z_steps_2 = 3;
d_dopant_to_CNT = 3.4;		% distance in y-direction [Ang]
rel_pos = 0;				% rel_pos = 0: dopant centered on honeycomb, rel_pos = 1: centered on C atom;

if strcmp(dopant_type,'KAuBr4')
	coords = KAuBr4_atoms_avg_relaxed_coord(file_to_read);
	dopant_molecule_1 = single_KAuBr4_gen(coords,rel_pos,sheet,[x_steps_1;z_steps_1],d_dopant_to_CNT);
	dopant_molecule_2 = single_KAuBr4_gen(coords,rel_pos,sheet,[x_steps_2;z_steps_2],d_dopant_to_CNT);
	
elseif strcmp(dopant_type,'AuBr4')
	coords = AuBr4_atoms_avg_relaxed_coord(file_to_read);
	dopant_molecule_1 = single_AuBr4_gen(coords,rel_pos,sheet,[x_steps_1;z_steps_1],d_dopant_to_CNT);
	dopant_molecule_2 = single_AuBr4_gen(coords,rel_pos,sheet,[x_steps_2;z_steps_2],d_dopant_to_CNT);
	
else	% K
	coords = K_atoms_avg_relaxed_coord(file_to_read);
	dopant_molecule_1 = single_K_gen(coords,rel_pos,sheet,[x_steps_1;z_steps_1],d_dopant_to_CNT);
% 	dopant_molecule_2 = single_K_gen(coords,rel_pos,sheet,[x_steps_2;z_steps_2],d_dopant_to_CNT);
end

dopant_coords_1 = num2cell(dopant_molecule_1.atoms);
dopant_coords_2 = num2cell(dopant_molecule_2.atoms);

%% Save to file and open in xcrysden
[hour,min,sec] = hms(datetime('now'));
[year,mon,day] = ymd(datetime('now'));

cd /home/khaiyi/Documents/Data/CNT_with_KAuBr4/

hour = num2str(hour);
min = num2str(min);
sec = num2str(round(sec));
year = num2str(year);
mon = num2str(mon);
day = num2str(day);
chirality = [num2str(chirality(1)),num2str(chirality(2))];

folder_name = [year,mon,day];
filename = [hour,':',min,':',sec,'_',chirality,'unrolled'];
if ~exist(['coord/',folder_name], 'dir')
  mkdir(['coord/',folder_name]);
end

% generate a 'VMD' coordinate file
outxyz([folder_name,'/',filename],sheet,dopant_coords_1{:},dopant_coords_2{:});
cmd = sprintf('xcrysden --xyz coord/%s/%s.xyz &',folder_name,filename);
system(cmd);