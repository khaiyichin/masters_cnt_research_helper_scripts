clear
close all
%% Define user parameters
directory = '/home/khaiyichin/masters_cnt_research';
cd(directory);
filename = '80SWNT_K_1_1x1x1_kgrid';
extension = '.AVTRANS';
file = [filename,extension];

interpolate = 0;	% perform linear interpolation? No:0; Yes:1
mplot = 1;			% plot graphs? No:0; Yes:1
mlabel = 1;			% display label on zero point? No:0; Yes:1
msave = 1;			% save all the figures? No:0; Yes:1
custom_title = 1;	% customize plot title? No:0; Yes:1
mdataout = 0;			% save zero point? No:0; Yes:1
stampede2 = 0;		% done on stampede2? No:0; Yes:1

x_trans = 0/10;
y_trans = 1.5/10;

if custom_title == 1
	plot_title = '(8,0) K-doped SWNT';
else
	plot_title = filename;
end

if stampede2 == 1
	factor = 2;
else
	factor = 1;
end

%% 
if  mdataout == 1
    fid_data = fopen('zero_data.txt','a');
    fprintf(fid_data,...
        '#\t%s\n case_name\t\tE /[eV]\t\tTrans /[G0]\n',...
        datetime);
end 


column_array_data = importdata(file);

if isa(column_array_data,'struct') == true		% when importing data, .TRANS gives a struct while .AVTRANS gives a double array
	x = column_array_data.data(:,1);
	y = column_array_data.data(:,2);
else
	x = column_array_data(:,1);																																																																																																																																																																																																																																																																																																																																																																																																																																																																																				
	y = column_array_data(:,2);
end

ind = find(x == 0);
ind = round(length(x)/2);
zero_x = x(ind);
zero_y = y(ind);

x_points = x;
y_points = y;

if interpolate == 1
	new_interval = abs((x(1)-x(2))/2);
	
	query_pts_x = [min(x):new_interval:0,0:new_interval:max(x)]';	% the interpolated points which you want to sample at (needs to be smaller than the current interval)
	query_pts_y = interp1(x,y,query_pts_x);
	
	ind = find(query_pts_x(:) == 0);
	zero_x = query_pts_x(ind);
	zero_y = query_pts_y(ind);
	
	x_points = query_pts_x;
	y_points = query_pts_y;
end

if mplot==1
	figure;
	plot(x_points,y_points/factor,'LineWidth',1.5);
	hold on
	grid on
	xlabel('E - E_F (eV)' ,'FontSize',15);
	ylabel('G/G_0','FontSize',15);
	set(gca,'FontSize',15);
	xticks(-2:0.5:2)
	yticks(0:1:9.5)
	axis([-2,2,0,9.5])
% 	yticks(0:0.5:2)
% 	axis([-2,2,0,2])
	if mlabel==1
		text(zero_x+x_trans,zero_y/factor+y_trans,sprintf('(%.3f, %.3f)',0,zero_y/factor),'FontSize',16);
		plot(0,zero_y/factor,'*r')
	end

	title(plot_title,'FontSize',14);%,'Interpreter','none');
	fig = gcf;
	fig.PaperUnits = 'inches';
	fig.PaperPosition = [0 0 6.75 5.5];

	if msave==1
		print(filename,'-dpng','-r0');
	end
end

if  mdataout==1
	fprintf(fid_data,'%s\t\t0\t\t%.7f\n',plot_title,zero_y);
end
