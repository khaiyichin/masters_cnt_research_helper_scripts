function rotation_matrix = dopant_orientation(rot_ang_x,rot_ang_y,rot_ang_z)
%dopant_orientation Generates rotation matrix
%    Rotates in the order of: about the x-axis, about the y-axis, about the
%    z-axis. Note that the CNT axis is along the z-axis. All arguments in
%    degrees

rot_ang_x = deg2rad(rot_ang_x);
rot_x = [1 0 0; 0 cos(rot_ang_x) -sin(rot_ang_x); 0 sin(rot_ang_x) cos(rot_ang_x)];
rot_ang_y = deg2rad(rot_ang_y);
rot_y = [cos(rot_ang_y) 0 sin(rot_ang_y); 0 1 0; -sin(rot_ang_y) 0 cos(rot_ang_y)];
rot_ang_z = deg2rad(rot_ang_z);
rot_z = [cos(rot_ang_z) -sin(rot_ang_z) 0; sin(rot_ang_z) cos(rot_ang_z) 0; 0 0 1];

rotation_matrix = rot_z*rot_y*rot_x;

end

