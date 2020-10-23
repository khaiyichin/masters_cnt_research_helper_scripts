function coords = ref_atom_coords_gen(chirality,box,CNT,position,num_dopants,d_dopants,offset,dopant_dist)
%ref_atom_coords_gen Generate coordinates for center atom of dopant
%   Strictly for armchair or zig-zag configuration only; for other 
%	configuration these coordinates	might be incorrect

xBox = box(1);
yBox = box(2);
d_CC = 1.421;

if position == 1
	sign = +1;
elseif position == -1
	sign = -1;
elseif position == 0
	sign = 0;
	dopant_dist = 0;
end

if chirality(1) == chirality(2) % armchair config
	x_coord = xBox/2;
	y_coord = yBox/2 + sign*(CNT.diameter/2 + dopant_dist);
	z_coord = sqrt(3)/2*d_CC;
	
	coords = [];
	for i = offset:d_dopants:(num_dopants*d_dopants)-(d_dopants-offset)
		coords = [	coords;
					x_coord, y_coord, (2*i-1)*z_coord	];
	end	
else % zig-zag config
	x_coord = xBox/2;
	y_coord = yBox/2 + sign*(CNT.diameter/2 + dopant_dist);
	z_coord = d_CC;
	
	coords = [];
	for i = offset:d_dopants:(num_dopants*d_dopants)-(d_dopants-offset)
		coords = [	coords;
					x_coord, y_coord, (3*i-2)*z_coord	];
	end	
end

end

