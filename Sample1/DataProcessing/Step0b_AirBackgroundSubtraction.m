% Use this script to subtract air background patterns from raw 2D
% diffraction patterns, scaled by the relative intensity of the patterns.

% These analyses support the publication of:
% "Mesoscale Structure and Composition Varies Systematically in Human Tooth
% Enamel", by R. Free, K. DeRocher, V. Cooley, R. Xu, S.R. Stock, and D. Joester.
% This script also includes the units and axes information for each plot.

% Author: Robert Free
% Editor: Derk Joester

%% Clean up
clear all;
close all;
clc
%% Flags
save_maps = false;
save_results = false;
plot_maps = true;
plot_WH = false;
%% Initialize Structs to hold data
Pattern=struct;
AzPlots=struct;
%% Define Parameters

%define detector axes
xAxis = linspace(1,2048,2048);
yAxis = linspace(2048,1,2048);
[Y,X] = ndgrid(yAxis,xAxis);
inner_lim = 150; %inner pixel limit for transform plot relative to beam center
outer_lim = 850; %outer pixel limit for transform plot relative to beam center

%define aspect ratio for plots and commonly used items for plots
xAsp=1;
yAsp=1;
XlabelStr='X-position (0.5um steps)';
YlabelStr='Y-position (0.5um steps)';
figWidth=500;
figHeight=300;

%Define dimensions of map in pixels (different for each scan)
xPts=52;
yPts=15;

%number of patterns
numPatt=xPts*yPts;
pattRange=1:numPatt;
pattMatrix=reshape(pattRange,[xPts,yPts])';
% pattSampleMatrix=pattMatrix(:,:); %Cropped area based on examining results
% pattSample=reshape(pattSampleMatrix',[1,numel(pattSampleMatrix)]); %Create list of patterns determined to be sample
% pattPlatinum=sort(pattMatrix(PtPatterns))';
% pattPlatinum=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,38,39,40,41,42,43,44,45,46,47,48,49,50,51,75,76,77,78,79,80,81,82,83,84,85,112,113,114,115,149,150,151,186,187,778,815,852,889,890];

num2Fit=4; %reflections in quadruplet

%experimental parameters
energyKeV=17.00000; %energy for scan
h=4.135667*10^-18; %keV*s
c=2.998792*10^18; %angstrom/s
lambda=h*c/energyKeV; %anstroms

%detector parameters
samp2det = 155.3818; %sample to detector distance in mm (may change per scan)
pixelSize = 0.079278; %pixel size in mm (usually the same for MAR detector from 34ide)
Sdd = samp2det / pixelSize; %sample to detector distance in pixels
tilt.rot = 167; % rotation of tilt plane relative to x-axis in degrees
tilt.angle = 0.675; % rotation of detector normal within the tilt plane in degrees
% xc = 1050.367; %refined fit of beam center in x direction from ceria (indexed from left side of pattern)
% yc = 888.856; %refined fit of beam center in y direction from ceria (indexed from top of pattern)
xc = 1051.623; %refined fit from pattern number 15, 121 reflection.
yc = 888.991; %refined fit from pattern number 15, 121 reflection.

DX = X-xc;
DY = Y-yc;

NtwoTheta = 1482; %Define the number of 2-theta bins
Nwedges = 120; %Define an EVEN number of azimuthal wedges for integration
window = 9; %Define an ODD integer for the number of azimuthal bins to span each wedge
Nchi = Nwedges*window;
twoThetaAxis = linspace(atand(inner_lim/Sdd),atand(outer_lim/Sdd),NtwoTheta);
chiAxis = linspace(0,360,Nchi+1); 
chiAxis = chiAxis(1:end-1);
Ntwotheta = numel(twoThetaAxis);

Kfactors=[1,1,1,1]; %Define Scherrer factors for each reflection

%functions relating angular width on detector at different radial positions
delta2theta = @(radDis) atand((radDis+(pixelSize/2))/samp2det)-atand((radDis-(pixelSize/2))/samp2det);
radialPixelDistance = @(twoTheta) Sdd*tand(twoTheta);

% define name of sample and scan (important for reading data)
sampleName='G2_OE3, ';
scanName='G2_OE3_fine';

%% Define regions for intensity normalization

Min2thetaBounds = 5.5;
Max2thetaBounds = 7.5;

fit2thetaMinInd=find(twoThetaAxis>5.5,1); %define index for left edge of quadruplet region (used for background subtraction)
fit2thetaMaxInd=find(twoThetaAxis<7.5,1,'last'); %and for right edge

%% Air Pattern Integration
%Read in pattern 1 (designated air pattern)
    filename=strcat(strcat(scanName,'_',num2str(1,'%03.0f'),'.tif'));
    AirRawImage=double(imread(filename));    
    %create interpolant from image based on best guess for beam center
    F = griddedInterpolant(flipud(DY),DX,flipud(AirRawImage),'linear');
    %Convert from polar to cartesian coordinates (2theta)
    queryPts = diffpol2cart_getQueryPts(twoThetaAxis,chiAxis,[0,0],[tilt.rot,tilt.angle],Sdd);
    AirCart2theta = reshape(F(queryPts(2,:),queryPts(1,:)),Nchi,NtwoTheta);
    %Compute 1D diffraction profile vs 2-theta
    AirDiffProfile=sum(AirCart2theta,1);
    %Integrate over designated background range
    AirIntensityMetric=sum(AirDiffProfile(fit2thetaMinInd:fit2thetaMaxInd));

%% Air Pattern Subtraction
IntMatrix=zeros(yPts,xPts);
ScaleMatrix=zeros(yPts,xPts);
%Loop over patterns in sample region
for j=pattRange
    %read in background corrected diffractin pattern .tif
    filename=strcat(strcat(scanName,'_',num2str(j,'%03.0f'),'.tif'));
    rawImage=double(imread(filename));    
    
    %create interpolant from image based on best guess for beam center
    F = griddedInterpolant(flipud(DY),DX,flipud(rawImage),'linear');
    
    %Convert from polar to cartesian coordinates (2theta)
    queryPts = diffpol2cart_getQueryPts(twoThetaAxis,chiAxis,[0,0],[tilt.rot,tilt.angle],Sdd);
    cart2theta = reshape(F(queryPts(2,:),queryPts(1,:)),Nchi,NtwoTheta);
    
    %Compute average intensity (2theta) in each pattern
    Pattern(j).AvgInt2theta=sum(sum(cart2theta,1),2)/(NtwoTheta*Nchi);   
    
    %Compute 1D diffraction profile vs 2-theta
    Pattern(j).DiffProfile=sum(cart2theta,1);
    
    %Integrate background range
    Pattern(j).IntensityMetric=sum(Pattern(j).DiffProfile(fit2thetaMinInd:fit2thetaMaxInd));
    Pattern(j).IntensityScale=Pattern(j).IntensityMetric/AirIntensityMetric;

    %Scale Air Pattern to current pattern and subtract
    ScaledAirImage=Pattern(j).IntensityScale*AirRawImage;
    BGcorrectedImage=uint16(rawImage-ScaledAirImage);

    %Write BG corrected Image to tiff
    filenameOut=strcat(strcat('scaled_corrected/',scanName,'_',num2str(j,'%04.0f'),'.tif'));
    imwrite(BGcorrectedImage,filenameOut,"Compression","none");
    
    %Assemble Matrices to visualize across all patterns for each beam
    %center under analysis
    IntMatrix(ceil(j/xPts),(1+mod(j-1,xPts)))=Pattern(j).AvgInt2theta;
    ScaleMatrix(ceil(j/xPts),(1+mod(j-1,xPts)))=Pattern(j).IntensityScale;
   
    disp(j)
end

%% Plot maps
if plot_maps
    fig=figure('Name','Int Matrix');
    set(fig,'Position',[200,200,500,300])
    imagesc(IntMatrix(:,:))
    daspect([yAsp,xAsp,1])
    colormap(gca,parula);
    title('Intensity Matrix')
    xlabel(XlabelStr)
    ylabel(YlabelStr)
    cbar = colorbar;
    caxis([0,1000])
    set(gca,'TickDir','out');
    cbar.TickDirection = 'out';

    fig=figure('Name','Scale Matrix');
    set(fig,'Position',[200,200,500,300])
    imagesc(ScaleMatrix(:,:))
    daspect([yAsp,xAsp,1])
    colormap(gca,parula);
    title('Scale Matrix')
    xlabel(XlabelStr)
    ylabel(YlabelStr)
    cbar = colorbar;
    caxis([0,1.1])
    set(gca,'TickDir','out');
    cbar.TickDirection = 'out';
end