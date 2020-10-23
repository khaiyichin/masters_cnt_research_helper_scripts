classdef molecule < handle
	%MOLECULE Dopant molecule which contains all 'dopant' objects (atoms)
	%   Inherits keeps track of all the 'dopant' atoms and coordinates
	
	properties
		atoms;					% All atoms in the dopant molecule; to be passed in as array; center_atom must be the first
		center_atom;			% The atom used as the reference point for other atoms in the molecule
		ref_coords;				% Coordinates of primary configuration of atoms in one molecule; must be in order as 'atoms_array' as a nx3 vector
	end
	
	methods
		% Constructor
		function obj = molecule(atoms_array,center_atom)
			obj.atoms = atoms_array;
			obj.center_atom = center_atom;
		end
		
		function relative_pos(obj,coords)	% 'coords' is nx3
			obj.ref_coords = coords;
		end
		
		function center_atom_duplication(obj,dupl_coord_array) % 'dupl_coord_array' is mx3
			obj.center_atom.coord = [obj.center_atom.coord;dupl_coord_array];
		end
		
		function rem_atoms_duplication(obj)
			dupl_coords = zeros(length(obj.atoms),3,length(obj.center_atom.coord(:,1)));
			
			for i = 1:length(obj.center_atom.coord(:,1))
				x = obj.center_atom.coord(i,1);
				y = obj.center_atom.coord(i,2);
				z = obj.center_atom.coord(i,3);
				
				dupl_coords(:,1,i) = obj.ref_coords(:,1) + x;
				dupl_coords(:,2,i) = obj.ref_coords(:,2) + y;
				dupl_coords(:,3,i) = obj.ref_coords(:,3) + z;
				
			end
			
			for i = 2:length(obj.atoms)
				obj.atoms(i).coord = reshape(dupl_coords(i,:,:),...
					[3,length(obj.center_atom.coord(:,1))])';	% reshaped coordinates is mx3
			end
		end
	end
	
end

