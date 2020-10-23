function obj = xyMove(obj,x,y)

obj.coord(:,1) = obj.coord(:,1) + x;
obj.coord(:,2) = obj.coord(:,2) + y;
