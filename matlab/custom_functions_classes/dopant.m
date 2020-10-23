classdef dopant < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        symbol
        coord
        unit
    end
    properties (Dependent)
        atomNo;
    end

    methods
        function value = get.atomNo(obj)
            value = size(obj.coord,1);
        end
    end

end
