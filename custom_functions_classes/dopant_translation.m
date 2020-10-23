function translation_matrix = dopant_translation(trans_x,trans_y,trans_z,coord_array_size)
%dopant_translation Generates coordinate translation matrix

translation_matrix = ...
	[	trans_x*ones(coord_array_size(1),1),...
		trans_y*ones(coord_array_size(1),1),...
		trans_z*ones(coord_array_size(1),1)	];

end

