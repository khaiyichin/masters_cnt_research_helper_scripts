function outxyz(varargin)
% the input sequence is (sys_name,molecule1,molecule1,molecule3 ....)


if ~exist('coord','dir')
    mkdir('coord');
end


if nargin <= 1
    msg ='Not enough input. The input sequence is (outputName,mole1,mole2...)';
    error(msg)
else
    n = nargin;
    file = varargin{1};
    addr = ['coord/',file,'.xyz'];

    total = 0;
    for i = 2:n
        row(i-1) = size(varargin{i}.coord,1);
        total = total + row(i-1);
    end

    % output the xyz file
    fileID = fopen(addr,'w');
    fprintf(fileID,' %d\n\n',total);

    for j=2:n
      for i = 1:row(j-1)
        fprintf(fileID,'%s\t%.8f\t%.8f\t%.8f \n',...
        varargin{j}.symbol, ...
        varargin{j}.coord(i,1), ...
        varargin{j}.coord(i,2), ...
        varargin{j}.coord(i,3));
      end
    end

    fclose(fileID);
end
