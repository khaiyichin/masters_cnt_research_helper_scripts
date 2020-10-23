function obj = zrotate(obj,pt1,pt2)
% obj = zrotate(obj,pt1,pt2)

ab = pt1 - pt2;

vector_x = [1 0];
theta = -acos(ab*vector_x'/norm(ab));
R = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];

obj.coord = (R*obj.coord')';
