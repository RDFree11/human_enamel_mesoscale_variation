% Use this script in conjunction with the .csv data files to generate rough
% versions of the figures for 
% "Mesoscale Structure and Composition Varies Systematically in Human Tooth
% Enamel", by R. Free, K. DeRocher, V. Cooley, R. Xu, S.R. Stock, and D. Joester.
% This script also includes the units and axes information for each plot.

% Author: Robert Free
% Editor: Derk Joester

%% Clean Up
clear variables 
close all
clc

%% Flags
save_figure = false;

%% Define units and figure parameters
% define resolution for d-spacing and azimuthal
dRes        = 2048; 
chiRes      = 2048;
twoThetaRes = 2048;

% define aspect ratio for plots and commonly used items for plots
xAsp = 1;
yAsp = 1;

% define pixel pitch for CLAA map
dx = 0.5; % [µm]
dy = 0.5; % [µm]

% crop CLAA map to indices
x_range = 8:22;

% labels for CLAA map
XlabelStr='x [µm]';
YlabelStr='y [µm]';

% Define dimensions of map in pixels (different for each scan)
xPts = 23;
yPts = 41;

energyKeV = 17.00000; %energy for scan

% physical constants
h = 4.135667*10^-18;    % [keV*s]
c = 2.998792*10^18;     % [angstrom/s]
lambda = h*c/energyKeV; % [angstrom]

samp2det  = 109.2976; % sample to detector distance in mm (may change per scan)
pixelSize = 0.079278; % pixel size in mm (MAR165 detector from APS 34IDE)

% anonymous functions
delta2theta = @(radDis) atand((radDis+(pixelSize/2))/samp2det)-atand((radDis-(pixelSize/2))/samp2det);
radialPixelDistance = @(twoTheta) samp2det*tand(twoTheta);

% define bounds of d-spacing for cartesian transform axis (extract from 1D
% diffraction profile constructed using same pixel values as cartesian
% transform in each case)
dMin = 1.3367920; % minimum d-spacing in angstroms
dMax = 3.9993496; % maximum d-spacing in angstroms
twoThetaMin = 7.7357362E-03; % minimum 2theta in degrees
twoThetaMax = 3.1677841E+01; % maximum 2theta in degrees

% Define x and y axis of plots based on how they were constructed in fit2D
dStep        = (dMax-dMin)/(dRes-1);
twoThetaStep = (twoThetaMax-twoThetaMin)/(twoThetaRes-1);
chiStep      = 360/chiRes;
dAxis        = [dMin:dStep:dMax];
twoThetaAxis = [twoThetaMin:twoThetaStep:twoThetaMax];
chiAxis      = [0:chiStep:360-chiStep];

% Define number of reflections and patterns to analyze
num2Fit       = 4;
numPatt       = 943;
fullPattRange = [1:numPatt];
pattRange     = [1:numPatt];
PtEdgeList    = [23:23:943];

% Define name of sample and scan (important for reading data)
sampleName = 'Enamel D, ';
scanName   = '50s_pattern';

% Define plot options
fig_width     = 6.5;
fig_asp_ratio = 2.5;
fig_zoom      = 3;
fig_pos       = [1,8,1+fig_width*fig_zoom,fig_width*fig_zoom/fig_asp_ratio]; %[in]
offset        = 10000;

FontSize_Panel = 10*fig_zoom;
FontSize_Label = 7*fig_zoom;
FontSize_Tick  = 6*fig_zoom;

%% FIGURE 2

mfile_name = 'Fig2_plot.m';

% paths & filenames
pn_csv     = './'
fn_inp     = {'CLAA_Map',...
          'Rod_Head_Cartesian_Diffraction_Pattern',...
          'Interrod_Cartesian_Diffraction_Pattern',...
          'Rod_Head_1D_Diffraction_Plot',...
          'Rod_Head_1D_Diffraction_Plot',...
          'Rod_Head_Quadruplet_Fit_Data',...
          'Interrod_Quadruplet_Fit_Data'};

pn_out = './';
fn_out = 'Fig2.eps';
          
% Import data
for ii = 1:length(fn_inp)
    M{ii} = readmatrix([pn_csv,fn_inp{ii},'.csv']);
end

%% Plot
set(groot,'defaultLineLineWidth',0.5);

fig2 = figure('Units','inches','Position',fig_pos);

% plot CLAA matrix
ax1 = subplot(2,4,[1,5]);
imagesc([0,diff(x_range([1,end]))*dx],[0,(size(M{1},1)-1)*dy],M{1}(:,x_range));
text(0.025,0.96,'A','Units','normalized','FontSize',FontSize_Panel,'Color','w');
daspect([yAsp,xAsp,1])
colormap bone;
ax1.FontSize = FontSize_Tick;
xlabel(XlabelStr,'FontSize',FontSize_Label);
ylabel(YlabelStr,'FontSize',FontSize_Label);
caxis([0,8]);
ax1.TickDir = 'out';

ax2 = subplot(2,4,2)
imagesc(dAxis,chiAxis,M{2});
text(0.01,0.92,'B','Units','normalized','FontSize',FontSize_Panel,'Color','w')
colormap(ax2,turbo);
caxis([0,200])
axis square;
ax2.FontSize = FontSize_Tick;
ax2.TickDir = 'out';
xlabel('D-spacing [$\rm{\AA}$]','Interpreter','latex','FontSize',FontSize_Label)
ylabel('Azimuthal Angle $\chi$ [$^\circ$]','Interpreter','latex','FontSize',FontSize_Label)
cbar2 = colorbar;
cbar2.TickLabels = {};
cbar2.TickDirection = 'out';
cbar2.Label.String = 'Intensity [a.u]';
cbar2.Label.FontSize = FontSize_Label;

ax3 = subplot(2,4,3)
plot(dAxis,M{4});
text(0.01,0.92,'C','Units','normalized','FontSize',FontSize_Panel,'Color','k')
xlabel('D-spacing [$\rm{\AA}$]','Interpreter','latex','FontSize',FontSize_Label)
ylabel('Intensity [a.u.]','FontSize',FontSize_Label)
xlim([dAxis(1),dAxis(end)]) 
ax3.FontSize = FontSize_Tick;
ax3.TickDir = 'out';
ax3.YTickLabel = {};

ax4 = subplot(2,4,4)
xx=M{6}(:,1);
plot(xx,M{6}(:,2)+offset,'o','Color','Black')
text(0.01,0.92,'D','Units','normalized','FontSize',FontSize_Panel,'Color','k')
xlabel('D-spacing [$\rm{\AA}$]','Interpreter','latex','FontSize',FontSize_Label)
ylabel('Intensity [a.u.]','FontSize',FontSize_Label)
hold on
plot(xx,M{6}(:,3)+M{6}(:,4)+M{6}(:,5)+M{6}(:,6)+offset,'Color','Blue')
plot(xx,M{6}(:,3)+offset,'Color','#D95319')
plot(xx,M{6}(:,4)+offset,'Color','#EDB120')
plot(xx,M{6}(:,5)+offset,'Color','#7E2F8E')
plot(xx,M{6}(:,6)+offset,'Color','#77AC30')
stem(xx,M{6}(:,2)-(M{6}(:,3)+M{6}(:,4)+M{6}(:,5)+M{6}(:,6)),'.','Color','#4DBEEE')
legend('data','fit','121 fit','112 fit','030 fit','022 fit','residual','Location','northeastoutside','FontSize',5)
xlim([xx(1),xx(end)])
ylim([-10000,140000])
ax4.TickDir = 'out';
ax4.FontSize = FontSize_Tick;
ax4.YTickLabel = {};
hold off

ax5 = subplot(2,4,6)
imagesc(dAxis,chiAxis,M{3});
text(0.01,0.92,'E','Units','normalized','FontSize',FontSize_Panel,'Color','w')
xlabel('D-spacing [$\rm{\AA}$]','Interpreter','latex','FontSize',FontSize_Label)
ylabel('Azimuthal Angle $\chi$ [$^\circ$]','Interpreter','latex','FontSize',FontSize_Label)
axis square;
caxis([0,200])
colormap(gca,turbo);
ax5.TickDir = 'out';
ax5.FontSize = FontSize_Tick;
cbar5 = colorbar;
cbar5.TickLabels = {};
cbar5.TickDirection = 'out';
cbar5.Label.String = 'Intensity [a.u]';
cbar2.Label.FontSize = FontSize_Label;

ax6 = subplot(2,4,7)
plot(dAxis,M{5})
text(0.01,0.92,'F','Units','normalized','FontSize',FontSize_Panel,'Color','k')
xlabel('D-spacing [$\rm{\AA}$]','Interpreter','latex','FontSize',FontSize_Label)
ylabel('Intensity [a.u.]','FontSize',FontSize_Label)
xlim([dAxis(1),dAxis(end)])
ax6.TickDir = 'out';
ax6.FontSize = FontSize_Tick;
ax6.YTickLabel = {};

ax7 = subplot(2,4,8)
xx=M{7}(:,1);
plot(xx,M{7}(:,2)+offset,'o','Color','Black')
text(0.01,0.92,'G','Units','normalized','FontSize',FontSize_Panel,'Color','k')
xlabel('D-spacing [$\rm{\AA}$]','Interpreter','latex','FontSize',FontSize_Label)
ylabel('Intensity [a.u.]','FontSize',FontSize_Label)
hold on
plot(xx,M{7}(:,3)+M{7}(:,4)+M{7}(:,5)+M{7}(:,6)+offset,'Color','Blue')
plot(xx,M{7}(:,3)+offset,'Color','#D95319')
plot(xx,M{7}(:,4)+offset,'Color','#EDB120')
plot(xx,M{7}(:,5)+offset,'Color','#7E2F8E')
plot(xx,M{7}(:,6)+offset,'Color','#77AC30')
stem(xx,M{7}(:,2)-(M{7}(:,3)+M{7}(:,4)+M{7}(:,5)+M{7}(:,6)),'.','Color','#4DBEEE')
legend('data','fit','121 fit','112 fit','030 fit','022 fit','residual','Location','northeastoutside','FontSize',5)
xlim([xx(1),xx(end)])
ylim([-10000,140000])
ax7.TickDir = 'out';
ax7.FontSize = FontSize_Tick;
ax7.YTickLabel = {};
hold off

if save_figure
    saveas(gcf,[pn_out,fn_out],'epsc');
    disp(['Saved as: ',pn_out,fn_out])
else
    disp('Figure not saved!');
end

