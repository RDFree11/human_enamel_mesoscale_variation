% Use this script in conjunction with the .csv data files to generate rough
% versions of the figures and tables for 
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
mfile_name          = 'FigS4_plot.m';

%% Flags
save_figure = false;
save_tables = false;
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
x_range = 9:22;

% labels for CLAA map
XlabelStr='x [µm]';
YlabelStr='y [µm]';

% Define dimensions of map in pixels (different for each scan)
xPts = 23;
yPts = 41;

energyKeV = 17.00000; %energy for scan

%physical constants
h = 4.135667*10^-18; %keV*s
c = 2.99792*10^18; %angstrom/s
lambda = h*c/energyKeV; %angstroms

samp2det  = 109.2976; %sample to detector distance in mm (may change per scan)
pixelSize = 0.079278; %pixel size in mm (usually the same for MAR detector from 34ide)

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
dStep=(dMax-dMin)/(dRes-1);
twoThetaStep=(twoThetaMax-twoThetaMin)/(twoThetaRes-1);
chiStep=360/chiRes;
dAxis=dMin:dStep:dMax;
twoThetaAxis=twoThetaMin:twoThetaStep:twoThetaMax;
chiAxis=0:chiStep:360-chiStep;

% Define number of reflections and patterns to analyze
num2Fit=4;
numPatt=943;
fullPattRange=1:numPatt;
pattRange=1:numPatt;
PtEdgeList=[23:23:943];

% Define name of sample and scan (important for reading data)
sampleName = 'Enamel D, ';
scanName   = '50s_pattern';

% Define plot options
fig_width     = 6.5;
fig_asp_ratio = 2.5;
fig_zoom      = 2;
fig_pos = [1,8,1+fig_width*fig_zoom,fig_width*fig_zoom/fig_asp_ratio]; %[in]
offset  = 10000;

FontSize_Panel = 10*fig_zoom;
FontSize_Label = 7*fig_zoom;
FontSize_Tick  = 6*fig_zoom;

ang  = char(197); 
panl = {'A','B','C','D','E'};
varname = {'CLAA','\Delta_c/c','\Delta_a/a','\Deltax_Mg','\Deltax_CO3'};

yl ={4+[-2,2],29+[-2.5,2.5],9.447+[-0.0025,0.0025],6.880+[-0.015,0.015]};


label_str  = {'$\frac{\Delta CLAA}{\bar{CLAA}}\;[\%]$',...
              '$\frac{\Delta c}{\bar{c}}\;[\%]$',...
              '$\frac{\Delta a}{\bar{a}}\;[\%]$',...
              '\Deltax_{Mg} [at%]',...
              '\Deltax_{CO_3} [at%]'};
          
cbar_str ={{'',''},{'\DeltaX_{Mg}','[at%]'},{'\DeltaX_{CO_3}','[at%]'}};
alpha_val = 0.25;

%% Read in data
clear M
% paths & filenames
pn_csv = './figure source data/';
fn_inp = 'Fit_Quality_and_Derived_Quantities';

pn_out = './';
fn_out = {'FigS4.eps','FigS4_PValues.csv','FigS4_CorrelationCoefficients.csv'};

% Import data
ResultsMatrixRaw=readmatrix([pn_csv,fn_inp,'.csv']);

%Drop first column (pattern index)
ResultsMatrixClipped=ResultsMatrixRaw(:,3:end);

%Remove outliers and compute correlation
outlierMatrix=zeros(size(ResultsMatrixClipped));
ResultsMatrixClean=ResultsMatrixClipped;
for i=[1:4]
    [B,TF]=rmoutliers(ResultsMatrixClipped(:,i),'mean');
    outlierMatrix(:,i)=TF;
    ResultsMatrixClean(find(TF==1),:)=NaN;
end

%Plot full correlation comparison
figure(5);
[R,PValue,H]=corrplot(ResultsMatrixClean,'varNames',{'CLAA','121 size','a parameter','c parameter'});

if save_figure
    saveas(gcf,[pn_out,fn_out{1}],'epsc');
    disp(['Saved as: ',pn_out,fn_out])
else
    disp('Figure not saved!');
end

if save_tables
    writematrix(PValue,[pn_out,fn_out{2}])
    writematrix(R,[pn_out,fn_out{3}])
else
    disp('data not saved!')
end