% Use this script in conjunction with the .csv data files to generate rough
% versions of the figures for 
% "Mesoscale Structural Gradients in Human Tooth Enamel", by R. Free, 
% K. DeRocher, V. Cooley, R. Xu, S.R. Stock, and D. Joester.
% This script also includes the units and axes information for each plot.

% Author: Robert Free
% Editor: Derk Joester

%% Clean Up
clear variables 
close all
clc

%% Set path (run script rather than section) 
% set path
mfile_name          = 'FigS6_plot.m';

%% Flags
save_figure = false;

%% Metadata
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

% labels for CLAA map
XlabelStr = 'x [µm]';
YlabelStr = 'y [µm]';

% Define dimensions of map in pixels (different for each scan)
xPts = 23;
yPts = 41;

energyKeV = 17.00000; %energy for scan

%physical constants
h = 4.135667*10^-18; % keV*s
c = 2.99792*10^18;  % angstrom/s
lambda = h*c/energyKeV; % angstrom

samp2det  = 109.2976; % sample to detector distance in mm (may change per scan)
pixelSize = 0.079278; % pixel size in mm (usually the same for MAR detector from 34ide)

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
fullPattRange = 1:numPatt;
pattRange     = 1:numPatt;
PtEdgeList    = [23:23:943];

% Define name of sample and scan (important for reading data)
sampleName = 'Enamel D, ';
scanName   = '50s_pattern';

% Define plot options
fig_width     = 11;
fig_asp_ratio = 2;
fig_zoom      = 3;
fig_pos       = [1,8,1+fig_width*fig_zoom,fig_width*fig_zoom/fig_asp_ratio]; %[in]
offset        = 10000;

FontSize_Panel = 10*fig_zoom;
FontSize_Label = 7*fig_zoom;
FontSize_Tick  = 6*fig_zoom;

ang       = char(197); 
cmapstr   = {'bone','parula','parula','parula'};
panl      = {'A','B','C','D','E','F','G','H','I'};
clims     = {[0,8],[23,33],[9.438,9.448],[6.85,6.92]};
dlim_abs  = [[0,8];[25,33];[9.438,9.448];[6.85,6.92]];
cbar_str1 ={{'CLAA','[a.u.]'},{'s_{121}','[nm]'},{'a',['[',ang,']']},{'c',['[',ang,']']}};
cbar_str2 ={{'CLAA*',''},{'s_{121}*','[%]'},{'a*','[%]'},{'c*','[%]'}};
cbar_fmt  = {'','%2.1f','%1.2f','%1.2f'};
cbar_pos  = [0,0.5,0.01,0.3];
rnd       = [1,1,2,2];

%% paths & filenames
pn_csv = './';
fn_inp     = {'121size_BCvariance_1',...
              '121size_BCvariance_2',...
              '121size_BCvariance_3',...
              '121size_BCvariance_4',...
              '121size_BCvariance_5',...
              '121size_BCvariance_6',...
              '121size_BCvariance_7',...
              '121size_BCvariance_8',...
              '121size_BCvariance_9',...
              'aParameter_BCvariance_1',...
              'aParameter_BCvariance_2',...
              'aParameter_BCvariance_3',...
              'aParameter_BCvariance_4',...
              'aParameter_BCvariance_5',...
              'aParameter_BCvariance_6',...
              'aParameter_BCvariance_7',...
              'aParameter_BCvariance_8',...
              'aParameter_BCvariance_9',...
              'cParameter_BCvariance_1',...
              'cParameter_BCvariance_2',...
              'cParameter_BCvariance_3',...
              'cParameter_BCvariance_4',...
              'cParameter_BCvariance_5',...
              'cParameter_BCvariance_6',...
              'cParameter_BCvariance_7',...
              'cParameter_BCvariance_8',...
              'cParameter_BCvariance_9'};

pn_out = './';
fn_out = 'FigS6.eps';
          
% Import data
for ii = 1:27
    M(:,:,ii) = readmatrix([fn_inp{ii},'.csv']);
end

% select data and convert units
for i=1:9
    M(:,:,i) = M(:,:,i)/10; % convert from Angstrom to nm
end

%% Plot
close all

indexArray = [1,2,3;4,5,6;7,8,9];


hf_sub(1)=figure('Name','121 Size');
hp(1) = uipanel('Parent',hf_sub(1),'Position',[200 200 1000 1000]);
for ii=1:9
    subplot(3,3,ii,'Parent',hp(1))
    imagesc(M(:,:,indexArray(ii)))
    daspect([yAsp,xAsp,1])
    colormap(gca,parula);
    title(['121size_',num2str(indexArray(ii))],'interpreter','none')
%     xlabel(XlabelStr)
%     ylabel(YlabelStr)
    cbar = colorbar;
    caxis([23,33])
%     caxis('auto')
    set(gca,'TickDir','out');
%     cbar.TickDirection = 'out';
end
hf_sub(2)=figure('Name','a parameter');
hp(2) = uipanel('Parent',hf_sub(2),'Position',[200 200 1000 1000]);
for ii=1:9
    subplot(3,3,ii,'Parent',hp(2))
    imagesc(M(:,:,indexArray(ii)+9))
    daspect([yAsp,xAsp,1])
    colormap(gca,parula);
    title(['aParameter_',num2str(indexArray(ii))],'interpreter','none')
%     xlabel(XlabelStr)
%     ylabel(YlabelStr)
    cbar = colorbar;
    caxis([9.438,9.448])
%     caxis('auto')
    set(gca,'TickDir','out');
%     cbar.TickDirection = 'out';
end

hf_sub(3)=figure('Name','c parameter');
hp(3) = uipanel('Parent',hf_sub(3),'Position',[200 200 1000 1000]);
for ii=1:9
    subplot(3,3,ii,'Parent',hp(3))
    imagesc(M(:,:,indexArray(ii)+18))
    daspect([yAsp,xAsp,1])
    colormap(gca,parula);
    title(['cParameter_',num2str(indexArray(ii))],'interpreter','none')
%     xlabel(XlabelStr)
%     ylabel(YlabelStr)
    cbar = colorbar;
    caxis([6.85,6.92])
%     caxis('auto')
    set(gca,'TickDir','out');
%     cbar.TickDirection = 'out';
end

hf_main = figure(4);
set(hf_main,'Position',[200 200 1000 500])
npanels = numel(hp);
hp_sub = nan(1,npanels);

for idx = 1:npanels
    hp_sub(idx) = copyobj(hp(idx),hf_main);
    set(hp_sub(idx),'Position',[(idx-1)/npanels,0,1/npanels,1]);
end


close(hf_sub(1))
close(hf_sub(2))
close(hf_sub(3))
%%
if save_figure
    saveas(gcf,[pn_out,fn_out],'epsc');
    disp(['Saved as: ',pn_out,fn_out])
else
    disp('Figure not saved!');
end
