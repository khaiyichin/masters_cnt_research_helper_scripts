%% Calculates specific conductivity of CNT models in .xyz files
clear

folder_to_read = '/home/khaiyi/Documents/Data/CNT_with_KAuBr4/SWNT/8-0/80SWNT_K_05/TACC_calcs/zSmall/conduction_calcs/';
file_to_read = '80SWNT_K_05.xyz';

[species,x,y,z] = importRelaxedxyz([folder_to_read,file_to_read]);

angstrom = 1e-10;
dalton = 1.660539040e-27;	% 1 dalton (unified atomic mass units) in kgs from https://physics.nist.gov/cgi-bin/cuu/Value?tukg

% Atomic mass units from DOI 10.1515/pac-2015-0305
K_mass_da = 38.0893;		% mass of K atom in da
Au_mass_da = 196.97;		% mass of Au atom in da
Br_mass_da = 79.904;		% mass of Br atom in da
C_mass_da = 12.011;			% mass of K atom in da

total_mass_da = 0;
species_num = zeros(length(species),1);
for num_of_atoms = 1:length(species)
	switch species{num_of_atoms}
		case 'K'
			total_mass_da = total_mass_da + K_mass_da;
		case 'Au'
			total_mass_da = total_mass_da + Au_mass_da;
		case 'Br'
			total_mass_da = total_mass_da + Br_mass_da;
		case 'C'
			total_mass_da = total_mass_da + C_mass_da;
	end
end

disp('---------------------------')
disp(file_to_read)
disp('---------------------------')

total_mass_kg = total_mass_da*dalton;
disp('Total mass of CNT =')
disp([num2str(total_mass_kg),' kg,	',num2str(total_mass_da),' da'])
disp(' ')

total_length_A = max(z) - min(z);
total_length_m = total_length_A*angstrom;
disp('Total length of CNT =')
disp([num2str(total_length_m),' m,	',num2str(total_length_A),' A'])
disp(' ')

len_squared_over_mass = total_length_m^2/total_mass_kg;
disp('Length^2/Mass =')
disp([num2str(len_squared_over_mass),' m^2/kg']);
disp(' ')

conductance_quantum = 7.7480917310e-5;
scaled_cond = input(['Transmission value for ',file_to_read,' = ']);
spec_conductivity = scaled_cond*conductance_quantum*len_squared_over_mass;

disp('Specific conductivity =')
disp([num2str(spec_conductivity),' S m^2/kg'])