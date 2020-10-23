classdef unrolled_CNT
% The constructor argument is (width (circumference of CNT), length (tube
% length of CNT), center coords of graphene, type of CNT,)

    properties
		symbol = 'C';
        coord;
		chirality;
    end
    properties (Dependent)
        atomNo;
    end
    
    methods
        function obj = unrolled_CNT(nX,nZ,origin,type)
            % Total 4 input variables
            % nX  : width of graphene (number of unit cell in X direction)
            % nZ  : length of graphene (number of unit cell in Z direction)
            % xyzO: coordinate of the midpoint in the graphene
            % type: type of graphene ('armchair' or 'zigzag')

            % bondLength between C-C
            d_CC = 1.421;
			d_HC = 1.09;
            
            Cnumber=nX*nZ*4;
            
            %% build the unicell of two type of graphene
            switch type
                case 0
                    cellX = d_CC * 3^0.5;
                    cellZ = d_CC * 3;
                    
                    XOrigin = origin(1)-(nX * cellX -d_CC*3^0.5/2)/2;
                    YOrigin = origin(2);
                    ZOrigin = origin(3)-(nZ * cellZ)/2;
                    
                    base0 = [0  , 0, 0;
                        cellX/2, 0, d_CC/2;
                        cellX/2, 0, d_CC/2 + d_CC;
                            0  , 0, d_CC*2;];
                    base  = [base0(:,1)+XOrigin, base0(:,2)+YOrigin, base0(:,3)+ZOrigin, ];
                    
					obj.chirality = [8,0];
					
                case 1                   
                    cellX = d_CC * 3;
                    cellZ = d_CC * 3^0.5;
                    
                    XOrigin = origin(1)-(nX * cellX -d_CC)/2;
                    YOrigin = origin(2);
                    ZOrigin = origin(3)-(nZ * cellZ )/2;
                    
                    base0 = [0     , 0, 0;
                        d_CC/2      , 0, cellZ/2;
                        d_CC/2 + d_CC, 0, cellZ/2;
                        d_CC*2      , 0, 0];                    
                    base  = [base0(:,1)+XOrigin, base0(:,2)+YOrigin, base0(:,3)+ZOrigin, ];
					
					obj.chirality = [5,5];
					
                otherwise
                    error('Please input "armchair" or "zigzag"')
            end
                  
            %% generate the C atom coordinate
            xyzlist = [];
            for i = 0:nX-1
                for j = 0:nZ-1
                    new = [base(:,1) + i*cellX, base(:,2), base(:,3) + j*cellZ];
                    xyzlist = [xyzlist; new];
                end
            end

            %% generate the complete graphene coordinate
            obj.coord = coord_zsort(xyzlist);
        end
        
        
        %% get the total number of the atom in graphene
        function value = get.atomNo(obj)
            value = size(obj.coord,1);
        end
        
    end
end

