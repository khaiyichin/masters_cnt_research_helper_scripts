classdef nanotube
% The constructor argument is (symbol, chiral, bondLength, l).
% This class is used to built nanotube structure.
% Note: The code regrading the nanotube function comes from
% source code in vmd 'nanotube.tcl'.

  properties
    bondLength;
    chiral;
    symbol;
    % name;
    coord;
    length;
    diameter;
  end
  properties (Dependent)
    atomNo;

  end
  methods
    function value = get.atomNo(obj)
      value = size(obj.coord,1);
    end

    function obj = nanotube(symbol, chiral, bondLength, l)

    obj.bondLength = bondLength;
    obj.chiral = chiral;
    obj.symbol = symbol;
    obj.length = l;

    n = chiral(1);
    m = chiral(2);
%     l = length;
    a = bondLength;

    G = gcd(2*m+n,2*n+m);

    % chiral vector
    chiral = a * (3 * (m^2 + m*n + n^2))^0.5;
    r = chiral / 2 / pi;
    obj.diameter = 2*r;
    l_cell = 3^0.5 * chiral / G;

    n_cell = floor(l/l_cell);
    if n_cell == 0
        error('The length of CNT is smaller than one unit cell')
    end

    nchord = 2*(m^2 + m*n + n^2)^0.5;
    rfudge = pi/nchord/sin(pi/nchord);

    i_min = 0;
    i_max = (n + (n+2*m)/G);
    j_min = (-(2*n+m)/G);
    j_max = m;

    num = 1;

    for j=j_min:j_max
      for i = i_min:i_max
        xprime1 = 3*a^2*(i*(2*n+m)+j*(n+2*m)) / (2*chiral);
        yprime1 = 3*3^0.5*a^2*(i*m-j*n) / (2*chiral);

        xprime2 = xprime1 + 3*a^2*(m+n) / (2*chiral);
        yprime2 = yprime1 - a^2 * 3^0.5 * (n-m) / (2*chiral);

        phi1 = xprime1/r;
        phi2 = xprime2/r;

        x1 = rfudge*r*cos(phi1);
        x2 = rfudge*r*cos(phi2);
        y1 = rfudge*r*sin(phi1);
        y2 = rfudge*r*sin(phi2);
        z1 = yprime1;
        z2 = yprime2;

        if xprime1 >=0 & xprime1 < chiral & yprime1 >= 0 & yprime1 < l_cell
          base1(num,1) = x1;
          base1(num,2) = y1;
          base1(num,3) = z1;

          base2(num,1) = x2;
          base2(num,2) = y2;
          base2(num,3) = z2;

          num = num + 1;
        end
      end
    end

    xyzlist = [base1(:,1) base1(:,2) base1(:,3)];
    xyzlist = [xyzlist; base2(:,1) base2(:,2) base2(:,3)];

    for i = 1:n_cell-1
      xyzlist = [xyzlist ; base1(:,1) base1(:,2) base1(:,3)+i*l_cell];
      xyzlist = [xyzlist ; base2(:,1) base2(:,2) base2(:,3)+i*l_cell];
    end

    repeat = [];
    m = 1;
    while m ~= length(xyzlist)
      for n=m+1:length(xyzlist)
        flag = abs(xyzlist(m,:) - xyzlist(n,:));
        if flag(1)<1e-3 & flag(2)<1e-3 & flag(3)<1e-3
          repeat = [repeat n];
        end
      end
      m = m+1;
    end
    xyzlist(repeat,:) = [];

    z_min = min(xyzlist(:,3));
    xyzlist(:,3) = xyzlist(:,3) - z_min;

    obj.coord = coord_zsort(xyzlist);
  end
end
end
