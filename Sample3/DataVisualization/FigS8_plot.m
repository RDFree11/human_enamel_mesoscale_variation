% Use this script in conjunction with the .csv data files to generate rough
% versions of the figures and tables for 
% "Mesoscale Structure and Composition Varies Systematically in Human Tooth
% Enamel", by R. Free, K. DeRocher, V. Cooley, R. Xu, S.R. Stock, and D. Joester.
% This script also includes the units and axes information for each plot.

% Author: Derk Joester

% clean Up
clear variables 
close all
clc

% flags
save_fig = false;

% get path
mfile_name          = 'FigS8_plot.m';

% rel paths & filenames
pn_csv     = './figure source data/';
fn_inp     = 'Data_from_Daculsi_and_Kerebel_1978';

pn_out = './';
fn_out = 'FigS8.eps';

% import options
opts1 = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts1.DataLines = [2, 2];
DataLines2 = [3, 11];
opts1.Delimiter = ",";

% Specify column names and types
opts1.VariableNames = ["distance", "thickness", "width", "numberDensity"];
opts1.VariableTypes = ["char", "char", "char", "char"];
VariableTypes2      = ["double", "double", "double", "double"];

% Specify file level properties
opts1.ExtraColumnsRule = "ignore";
opts1.EmptyLineRule = "read";

opts2 = opts1;
opts2.DataLines = DataLines2;
opts2.VariableTypes = VariableTypes2;

% Import data
U = readmatrix([pn_csv,fn_inp,'.csv'],opts1);
T = readtable([pn_csv,fn_inp,'.csv'],opts2);


%% constants

% set half angle a crystallite vertex in basal plane
ang = 60; % [˚]

T.area = T.thickness.*T.width-T.thickness.^2/2/tand(ang); % nm^2
T.volumeFraction = T.numberDensity.*T.area*1e-6; % [a.u.] 

%% fit

% fit only data below threshhold distance
d_max = 50;

% select data to fit
inda = ~isnan(T.distance) & T.distance <= d_max;
% fit cross-sectional area vs. distance with linear model
fa = fit(T.distance(inda),T.area(inda),'poly1');
% prediction interval for fitted function, non-simultaneous
pa = predint(fa,T.distance(inda),0.95,'functional','off');
% confidence interval for fit parameters
CIa = confint(fa);

indb = ~isnan(T.distance) & ~isnan(T.volumeFraction);
% fit cross-sectional area vs. distance with linear model
fb = fit(T.distance(indb),T.volumeFraction(indb),'poly1');
% prediction interval for fitted function, non-simultaneous
pb = predint(fb,T.distance(indb),0.95,'functional','off');
% confidence interval for fit parameters
CIb = confint(fb);

% define plot axes
dav = [0,d_max];
dbv = T.distance([find(indb,1,'first'),find(indb,1,'last')]);

%% plot 
figure('Units',"inches","Position",[1,1,8,4]);
ax1 = subplot(1,2,1);


plot(dav,fa(dav),'r:'); hold on;
plot(T.distance(inda),T.area(inda),'o','MarkerEdgeColor','none','MarkerFaceColor','b'); hold on;
plot(T.distance(inda),pa,'r--','DisplayName','CI(95%)')

xlabel('Distance [µm]');
ylabel('Area [nm^2]');
legend('Location',"best")
ax1.TickDir = 'out';
text(d_max/2,fa(d_max),{'$\hat{A} = \hat{p}_1 \cdot d + \hat{p}_2$',...
               ['$\hat{p}_1 = ',num2str(fa.p1,'%2.2f'),'\;(',num2str(CIa(1,1),'%2.2f'),',',num2str(CIa(2,1),'%2.2f'),')','\;\rm{nm}  $'],...
               ['$\hat{p}_2 = ',num2str(fa.p2,'%2.2f'),'\;(',num2str(CIa(1,2),'%2.2f'),',',num2str(CIa(2,2),'%2.2f'),')','\;\rm{nm}^2 $']},...
               'Interpreter',"latex",'HorizontalAlignment',"center",'VerticalAlignment',"bottom");
grid on
ax1.YMinorTick = 'on';

ax2 = subplot(1,2,2);
plot(dbv,fb(dbv),'r:'); hold on;
plot(T.distance(indb),T.volumeFraction(indb),'o','MarkerEdgeColor','none','MarkerFaceColor','b'); hold on;
%plot(T.distance(indb),pb,'r--','DisplayName','CI(95%)');
text(250,0.9,{'$\hat{A} = \hat{p}_1 \cdot d + \hat{p}_2$',...
               ['$\hat{p}_1 = ',num2str(fb.p1,'%2.2g'),'\;(',num2str(CIb(1,1),'%2.2g'),',',num2str(CIb(2,1),'%2.2g'),')','\;\rm{nm}  $'],...
               ['$\hat{p}_2 = ',num2str(fb.p2,'%2.2g'),'\;(',num2str(CIb(1,2),'%2.2g'),',',num2str(CIb(2,2),'%2.2g'),')','\;\rm{nm}^2 $']},...
               'Interpreter',"latex",'HorizontalAlignment',"center",'VerticalAlignment',"bottom");

xlabel('Distance [µm]');
ylabel('Volume Fraction [a.u.]');
legend('Location',"southeast")
ax2.TickDir = 'out';
ax2.YMinorTick = 'on';
ax2.XMinorTick = 'on';
grid on
ylim([0,1.1])

if save_fig
    saveas(gcf,[pn_out,fn_out],'epsc');
    disp(['Saved as: ',pn_out,fn_out])
else
    disp('Figure not saved!');
end
