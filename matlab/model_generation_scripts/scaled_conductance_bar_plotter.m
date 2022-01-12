%% Scaled specific conductance calculator
% Generates a bar chart for a specific doped CNT (both semiconducting and
% metallic.

clear
close all

%%%%%%%%%%%%%%%%%%%%%%% USER CONFIG %%%%%%%%%%%%%%%%%%%%%%%
output_folder = '/home/khaiyichin/Downloads/cnt/';

filename_1 = 'k_doped_cond';
filename_2 = 'k_doped_scaled_cond';

ext_1 = '.fig';
ext_2 = '.png';

output_fig_1 = strcat(output_folder, filename_1, ext_1);
output_fig_2 = strcat(output_folder, filename_2, ext_1);
output_png_1 = strcat(output_folder, filename_1, ext_2);
output_png_2 = strcat(output_folder, filename_2, ext_2);

grayscale_image = 1; % 1 for grayscale, 0 for color

% Dopant mass
KAuBr4 = 38.0893 + 196.97 + 4*79.904;
K = 38.0893;
AuBr4 = 196.97 + 4*79.904;

% Transmission values G/G_0
trans_1 = 1.644;		% (8,0) lower dopant concentration
trans_2 = 2.958;		% (5,5)	lower dopant concentration
trans_3 = 2.181;		% (8,0) higher dopant concentration
trans_4 = 3.8385;		% (5,5) higher dopant concentration

% Number of dopants in one SWNT unit cell
duc_1 = 0.5;
duc_2 = 0.5;
duc_3 = 1;
duc_4 = 1;

% Mass of each dopant
m_d_1 = K;
m_d_2 = K;
m_d_3 = K;
m_d_4 = K;

%%%%%%%%%%%%%%%%%%%%%%% END OF USER CONFIG %%%%%%%%%%%%%%%%%%%%%%%

% Number of C atoms in the computational unit cell
C_uc_1 = 32;
C_uc_2 = 20;
C_uc_3 = 32;
C_uc_4 = 20;

% Mass of each C atom
m_C = 12.011;

% SWNT unit cell length in angstroms
len_uc_1 = 3*1.421;
len_uc_2 = sqrt(3)*1.421;
len_uc_3 = 3*1.421;
len_uc_4 = sqrt(3)*1.421;

% Scaled specific conductance
scaled_1 = trans_1 / ((m_d_1*duc_1+m_C*C_uc_1)/len_uc_1);
scaled_2 = trans_2 / ((m_d_2*duc_2+m_C*C_uc_2)/len_uc_2);
scaled_3 = trans_3 / ((m_d_3*duc_3+m_C*C_uc_3)/len_uc_3);
scaled_4 = trans_4 / ((m_d_4*duc_4+m_C*C_uc_4)/len_uc_4);

x_label = {'0.50','0.50','1.00','1.00'};
x_data = [1,2,4,5];

% Bar chart for conductance
figure(1)
hold on
y_data_1 = [trans_1;trans_2;trans_3;trans_4];
for i = 1:4
	if i == 1 || i == 3
        if grayscale_image == 1
            color_str = [179/255, 179/255, 179/255]; % grayscale colors
        else
            color_str = 'r';
        end
    else
        if grayscale_image == 1
            color_str = [77/255, 77/255, 77/255]; % grayscale colors
        else
            color_str = 'b';
        end
	end
	bar(x_data(i),y_data_1(i),'FaceColor', color_str);
	text(x_data(i)-0.285,y_data_1(i)+0.15,num2str(y_data_1(i),'%.3f'),'FontSize',12);
end

fig_1 = gcf;
fig_1.PaperUnits = 'inches';
fig_1.PaperPosition = [0 0 6.75 5.5];

axes_1 = gca;
axes_1.FontSize = 12;
axes_1.XTick = x_data;
axes_1.XTickLabel = x_label;
axes_1.XLabel.String = 'Dopant Concentration (atoms/unit cell)';
axes_1.XLabel.Position = axes_1.XLabel.Position + [0 -0.0375 0];
axes_1.YLim = [0,5];
axes_1.YTick = [0:0.5:axes_1.YLim(end)];
axes_1.YGrid = 'on';
axes_1.Title.String = 'Conductance (K Dopant)';
axes_1.YLabel.String = 'G/G_0';
axes_1.Box = 'on';
legend({'(8,0) SWNT','(5,5) SWNT'},'Location','Northwest','FontSize',12)

plot(axes_1.XLim,[2 2],'k--')
axes_1.Legend.String(3) = {'Undoped (5,5) SWNT'};
savefig(fig_1,output_fig_1)
print(output_png_1,'-dpng','-r0')

% Bar chart for scaled specific conductance
figure(2)
hold on
y_data_2 = [scaled_1;scaled_2;scaled_3;scaled_4];
for i = 1:4
	if i == 1 || i == 3
        if grayscale_image == 1
            color_str = [179/255, 179/255, 179/255]; % grayscale colors
        else
            color_str = 'r';
        end
    else
        if grayscale_image == 1
            color_str = [77/255, 77/255, 77/255]; % grayscale colors
        else
            color_str = 'b';
        end
	end
	bar(x_data(i),y_data_2(i),'FaceColor', color_str);
	text(x_data(i)-0.35,y_data_2(i)+0.0015,num2str(y_data_2(i),'%.4f'),'FontSize',12);
end

fig_2 = gcf;
fig_2.PaperUnits = 'inches';
fig_2.PaperPosition = [0 0 6.75 5.5];

axes_2 = gca;
axes_2.FontSize = 12;
axes_2.XTick = x_data;
axes_2.XTickLabel = x_label;
axes_2.XLabel.String = 'Dopant Concentration (atoms/unit cell)';
% axes_2.XLabel.Position = axes_2.XLabel.Position + [0 -0.00085 0];
axes_2.XLabel.Position = axes_2.XLabel.Position + [0 -0.0002 0];
axes_2.YLim = [0,0.04];
axes_2.YTick = [0:0.005:axes_1.YLim(end)];
axes_2.YGrid = 'on';
axes_2.Title.String = 'Specific Conductance (K Dopant)';
axes_2.YLabel.String = ['LG/mG_0 (',char(197),' amu^{-1})'];
axes_2.Box = 'on';
legend({'(8,0) SWNT','(5,5) SWNT'},'Location','Northwest','FontSize',12)

plot(axes_2.XLim,[.0204 .0204],'k--')
axes_2.Legend.String(3) = {'Undoped (5,5) SWNT'};
savefig(fig_2,output_fig_2)
print(output_png_2,'-dpng','-r0')
% filename = '';
% print(filename,'-dpng','-r0');
