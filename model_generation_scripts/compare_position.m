%% Compare positions of two .xyz files
clear

cd /home/khaiyi/Documents/Data/CNT_with_KAuBr4/SWNT/8-0/80SWNT_K_1/TACC_calcs/cell_size_test/

filename_1 = 'variable_cell/80SWNT_K_1_var_trans.xyz';
filename_2 = 'non_variable_cell/80SWNT_K_1_nonvar.xyz';

[~,x_1,y_1,z_1] = importRelaxedxyz(filename_1);
[~,x_2,y_2,z_2] = importRelaxedxyz(filename_2);

points_1 = [x_1,y_1,z_1];
points_2 = [x_2,y_2,z_2];

distance = [];
for i = 1:length(x_1)
	distance(i) = norm(points_1(i,:)-points_2(i,:),2);
end

disp('Average displacement')
disp(mean(distance))