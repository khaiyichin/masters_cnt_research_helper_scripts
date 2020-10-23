%% Translate coordinates in relaxed .xyz file into new .xyz file
clear

cd /home/khaiyi/Documents/Data/CNT_with_KAuBr4/

file_read = 'SWNT/8-0/80SWNT_K_1/TACC_calcs/cell_size_test/variable_cell/80SWNT_K_1_var.xyz';
file_write = 'SWNT/8-0/80SWNT_K_1/TACC_calcs/cell_size_test/variable_cell/80SWNT_K_1_var_trans.xyz';
x_trans = 0.4225575;
y_trans = 0.40137;
z_trans = -2.37e-4;

[species,x,y,z] = importRelaxedxyz(file_read);

species_new = string(species);
x_new = x - x_trans;
y_new = y - y_trans;
z_new = z - z_trans;

fileID = fopen(file_write,'w');
fprintf(fileID,'  %s\n\n',num2str(length(species)));
fprintf(fileID,'%s      %3.6f   %3.6f   %3.6f\n',[species_new,x_new,y_new,z_new]');
fclose(fileID);