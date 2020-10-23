function [diam,circ] = cnt_diam_n_circ(n,m)
%cnt_diam_n_circ Calculates diameter and circumference of the carbon
%nanotube
%   Given the chirality vectors n and m, this function computes the
%   diameter and circumference of the carbon nanotube.

circ = 1.2306*sqrt(n^2 + n*m + m^2);			% 1.2306 is the length of the unit vector a1 == a2
diam = circ/pi;
disp(['Diameter = ',num2str(diam)]);
disp(['Circumference = ',num2str(circ)]);

end

