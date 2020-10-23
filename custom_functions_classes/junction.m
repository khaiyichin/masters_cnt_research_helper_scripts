classdef junction
% the upper left and lower right will be chopped
    properties
        unit;
        intTube;
        overLap;
        sysL;
        rawCoord;
        coord;
        symbol;
        diameter;
        upperOffset;
        upperEnd;
        lowerEnd;
    end

    properties (Dependent)
        atomNo;
    end

    methods
        function obj = junction(raw)
            obj.rawCoord = raw;
        end

        function obj = MakeJunc(obj)
            tempUp = obj.rawCoord;
            tempUp(:,2) = tempUp(:,2) + (obj.intTube + obj.diameter) / 2;
            tempUp(:,3) = tempUp(:,3) - obj.upperOffset;

%             indx = find(tempUp(:,3) < ((obj.sysL+obj.unit/2+obj.upperOffset)/2 - 0.1 - obj.overLap/2) ); % only for transiesta cnt55 junction
            indx = find(tempUp(:,3) < (((obj.sysL)+obj.upperOffset)/2 - 0.1 - obj.overLap/2) ); % for others
            
            tempUp(indx,:) = [];
            obj.upperEnd = min(tempUp(:,3));

            tempDown = obj.rawCoord;
            tempDown(:,2) = tempDown(:,2) - (obj.intTube + obj.diameter) / 2;

%             indx = find(tempDown(:,3) > (((obj.sysL+obj.unit/2)+obj.upperOffset)/2 + 0.1 + obj.overLap/2)); % only for transiesta cnt55 junction
            indx = find(tempDown(:,3) > (((obj.sysL)+obj.upperOffset)/2 - 0.1 + obj.overLap/2)); % for others
            tempDown(indx,:) = [];

            upperMax = max(tempUp(:,3));
            indx = find(tempDown(:,3) > upperMax);
            tempDown(indx,:) = [];
            obj.lowerEnd = max(tempDown(:,3));

            obj.coord = coord_zsort([tempDown; tempUp]);
        end

        function value = get.atomNo(obj)
            value = size(obj.coord,1);
        end
    end
end
