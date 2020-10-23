function output_struct = dopant_selector(dopant_type_string)
%dopant_selector Selects reference coordinates for dopants
%   The dopant coordinates are obtained from the relaxed configuration of 
%	the individual molecules, with respect to a 'center' atom at the origin
if strcmp(dopant_type_string,'KAuBr4')
	ref_coords = ...				% KAuBr4 dopant
		[	0.0000	0.0000	0.0000	;
			2.1616	-2.1373	2.7727	;
			-0.0134	-2.4093	0.0000	;
			2.4146	0.0000	0.0000	;
			0.0134	2.4093	0.0000	;
			-2.4146	0.0000	0.0000	];
	
	dopant_gen_fcn = 1;
		
elseif strcmp(dopant_type_string,'AuBr4')
	ref_coords = ...				% AuBr4 dopant
		[	0.0000	0.0000	0.0000	;
			-0.0134	-2.4093	0.0000	;
			2.4146	0.0000	0.0000	;
			0.0134	2.4093	0.0000	;
			-2.4146	0.0000	0.0000	];
		
	dopant_gen_fcn = 2;
	
elseif strcmp(dopant_type_string,'K')
	ref_coords = ...				% K dopant
		[	0.000	0.000	0.000	];
	
	dopant_gen_fcn = 3;
end

output_struct = struct('coords',ref_coords,'fcn_num',dopant_gen_fcn);

end

