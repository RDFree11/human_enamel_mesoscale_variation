% These analyses support the publication of:
% "Mesoscale Structural Gradients in Human Tooth Enamel", by R. Free, 
% K. DeRocher, V. Cooley, R. Xu, S.R. Stock, and D. Joester.
% This script also includes the units and axes information for each plot.

% Author: Robert Free
% Editor: Derk Joester

%%
clear all;
close all;
clc

%% Initialize Strcuts to hold data
Centers=struct;
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
figWidth=300;
figHeight=500;

%experimental parameters
energyKeV=17.10000; %energy for scan
h=4.135667*10^-18; %keV*s
c=2.99792*10^18; %angstrom/s
lambda=h*c/energyKeV; %anstroms

%detector parameters
samp2det = 155.6166; %sample to detector distance in mm (may change per scan)
pixelSize = 0.079278; %pixel size in mm (usually the same for MAR detector from 34ide)
Sdd = samp2det / pixelSize; %sample to detector distance in pixels
tilt.rot = 167.8; % rotation of tilt plane relative to x-axis in degrees
tilt.angle = 0.718; % rotation of detector normal within the tilt plane in degrees

% xc = 1055.635; %refined fit of beam center in x direction from ceria (indexed from left side of pattern)
% yc = 888.119; %refined fit of beam center in y direction from ceria (indexed from top of pattern)
xc = 1055.775; %refined fit of beam center in x direction from 121 of full enamel sample
yc = 889.044; %refined fit of beam center in y direction from 121 of full enamel sample

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

%functions relating angular width on detector at different radial positions
delta2theta = @(radDis) atand((radDis+(pixelSize/2))/samp2det)-atand((radDis-(pixelSize/2))/samp2det);
radialPixelDistance = @(twoTheta) Sdd*tand(twoTheta);

% define name of sample and scan (important for reading data)
sampleName='Victoria Enamel, ';
scanName='ThinEnamel_fine';
%% Define fitting functions (check working directory for referenced functions)

%set up combination of 4 psuedo-voigt profiles as anonymous function
PseudoVoigt4 = @(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,x) PseudoVoigt(x,a1,a2,a3,a4)+PseudoVoigt(x,a5,a6,a7,a8)+PseudoVoigt(x,a9,a10,a11,a12)+PseudoVoigt(x,a13,a14,a15,a16);
PseudoVoigt2 = @(a1,a2,a3,a4,a5,a6,a7,a8,x) PseudoVoigt(x,a1,a2,a3,a4)+PseudoVoigt(x,a5,a6,a7,a8);
PseudoVoigt1 = @(a1,a2,a3,a4,x) PseudoVoigt(x,a1,a2,a3,a4);

%define combined full width half max function for psuedo-voigt profile
FWHMcalc = @(fg,fl) ((fg.^5)+2.69269.*(fg.^4).*(fl)+2.42843.*(fg.^3).*(fl.^2)+4.47163.*(fg.^2).*(fl.^3)+0.07842.*(fg).*(fl.^4)+(fl.^5)).^(1/5);

%define etta computation for psuedo voigt
nCalc = @(fg,fl) (1.36603*(fl./FWHMcalc(fg,fl))-0.47719*(fl./FWHMcalc(fg,fl)).^2+0.11116*(fl./FWHMcalc(fg,fl)).^3);
%% Calculate instrumental broadening (check working directory for referenced functions and data)

% Import 2theta pattern
tmpData=importdata('Si_2theta.txt'); %creates file name string based on loop
tmpData=tmpData.data; %extracts only data from imported file
tmpData=tmpData(2:end);
tmpData=reshape(tmpData,[2,length(tmpData)/2]);
tmpData=transpose(tmpData);
SiCalibAxis=tmpData(:,1); %writes 2theta axis vector
SiCalibInt=tmpData(:,2); %writes Intensity

%plot pattern
figure;
title('Si Standard Diffraction Profile')
plot(SiCalibAxis,SiCalibInt)
xlabel('Bragg Angle (2-theta)')
ylabel('Intensity')

%subtract background
fit2thetaMinInd=find(SiCalibAxis>20.4,1); %define index for left edge of quadruplet region (used for background subtraction)
fit2thetaMaxInd=find(SiCalibAxis<21.2,1,'last'); %and for right edge
xDiff=fit2thetaMaxInd-fit2thetaMinInd;
yDiff=SiCalibInt(fit2thetaMaxInd)-SiCalibInt(fit2thetaMinInd);
newXaxis=[1:length(SiCalibAxis)]';
newXaxis=newXaxis-fit2thetaMinInd;
yAdjustment=(yDiff/xDiff)*newXaxis+SiCalibInt(fit2thetaMinInd);
AdjustedProfile=SiCalibInt-yAdjustment;
negativeIndices=find(AdjustedProfile<0); %find indices of all negative values
AdjustedProfile(negativeIndices)=0; %set negative values to 0
SiCalibIntAdj=AdjustedProfile; % write to struct

%fit peak in 2theta space and show fit quality on plot
[SiFit,SiFitGOF]=fit(SiCalibAxis(fit2thetaMinInd:fit2thetaMaxInd),SiCalibIntAdj(fit2thetaMinInd:fit2thetaMaxInd),PseudoVoigt1,'StartPoint',[55,20.74,.1,.1],'Lower',[0.1,0,.000001,.000001])
SiFWHM=FWHMcalc(SiFit.a3,SiFit.a4);

figure;
plot(SiFit,'Blue')
hold on
plot(SiCalibAxis(fit2thetaMinInd:fit2thetaMaxInd)',SiCalibIntAdj(fit2thetaMinInd:fit2thetaMaxInd),'o','Color','Black')
plot(SiFit,SiCalibAxis(fit2thetaMinInd:fit2thetaMaxInd)',SiCalibIntAdj(fit2thetaMinInd:fit2thetaMaxInd),'residuals')
xlabel('Bragg Angle (2-theta)')
ylabel('Intensity')
title('Pseudo-Voigt Fit to Si Standard')
legend('fitted curve','data','residual','zero line');
hold off

xxx = -2:0.001:2;
FitMax=PseudoVoigt(0,SiFit.a1,0,SiFit.a3,SiFit.a4);
leftEdgeIndex = find(PseudoVoigt(xxx,SiFit.a1,0,SiFit.a3,SiFit.a4)>(FitMax/2),1);
rightEdgeIndex = find(PseudoVoigt(xxx,SiFit.a1,0,SiFit.a3,SiFit.a4)>(FitMax/2),1,'last');
SiFWHM2 = xxx(rightEdgeIndex)-xxx(leftEdgeIndex);
%% Define fit parameters (for 2theta fits)

%specify peaks to fit
num2Fit=4;
PeakRef={'121','112','030','022'}; %names of peaks for generating plot titles

%specify initial values for fitting parameters (fit stability is sensitive
%to these. Determined through itereative fitting and checking quality of
%fit)
TwoThetaGuesses=[14.86,15.08,15.37,15.93]; %two-theta guesses for quadruplet
IntGuess=20000; %intensity guess
FWHMgGuess=0.08; %FWHM 2theta of gaussian component guess
FWHMlGuess=0.08; %FWHM 2theta of lorentzian component guess
%combine all initial guesses into StartVals vector for fit function
StartVals=[IntGuess,TwoThetaGuesses(1),FWHMgGuess,FWHMlGuess,IntGuess,TwoThetaGuesses(2),FWHMgGuess,FWHMlGuess,IntGuess,TwoThetaGuesses(3),FWHMgGuess,FWHMlGuess,IntGuess,TwoThetaGuesses(4),FWHMgGuess,FWHMlGuess]; %initial guesses for parameters

% Define limits for parameters to explore during fitting
twoThetaVariation=0.2; % max size of permissible variation in 2theta from initial guess
FWHMMin=min(SiFit.a3,SiFit.a4); %use the minimum from the silicon fit to determine minimum FWHM allowed during fitting.
FWHMMax=0.15; % max FWHM Variation (in either gaussian or lorentzian component)from initial guess allowed
IntMin=0.0001;
IntMax=10000000;
LowerBounds=[IntMin,TwoThetaGuesses(1)-twoThetaVariation,FWHMMin,FWHMMin,IntMin,TwoThetaGuesses(2)-twoThetaVariation,FWHMMin,FWHMMin,IntMin,TwoThetaGuesses(3)-twoThetaVariation,FWHMMin,FWHMMin,IntMin,TwoThetaGuesses(4)-twoThetaVariation,FWHMMin,FWHMMin];
UpperBounds=[IntMax,TwoThetaGuesses(1)+twoThetaVariation,FWHMMax,FWHMMax,IntMax,TwoThetaGuesses(2)+twoThetaVariation,FWHMMax,FWHMMax,IntMax,TwoThetaGuesses(3)+twoThetaVariation,FWHMMax,FWHMMax,IntMax,TwoThetaGuesses(4)+twoThetaVariation,FWHMMax,FWHMMax];

fit2thetaMinInd=find(twoThetaAxis>14.0,1); %define index for left edge of quadruplet region (used for background subtraction)
fit2thetaMaxInd=find(twoThetaAxis<16.3,1,'last'); %and for right edge
%% Define beam center and beam center sampling matrices

Ncenter_sampling = 1;
center_range = 0.5;

xoffset_array = linspace(-center_range,center_range,Ncenter_sampling);
yoffset_array = linspace(-center_range,center_range,Ncenter_sampling);

xcenters_array = xoffset_array+xc;
ycenters_array = yoffset_array+yc;

[Xoffsets,Yoffsets]=meshgrid(xoffset_array,yoffset_array);
[Xcenters,Ycenters]=meshgrid(xcenters_array,ycenters_array);

DX = X-xc;
DY = Y-yc;

numCenters = numel(Xcenters);
%% Read in Selected Pattern
j=90; %Select index of file of interest

filename='AVG_EntireSample.tif';
% filename=strcat(strcat(scanName,'_',num2str(j,'%04.0f'),'.tif'));
rawImage=double(imread(filename));

figure;
imagesc(xAxis,yAxis,rawImage);
colormap('turbo')
caxis([0,150])
daspect([1,1,1])

F = griddedInterpolant(flipud(DY),DX,flipud(rawImage),'linear');
queryPts = diffpol2cart_getQueryPts(twoThetaAxis,chiAxis,[0,0],[tilt.rot,tilt.angle],Sdd);
cart2theta = reshape(F(queryPts(2,:),queryPts(1,:)),Nchi,NtwoTheta); 
DiffProfile=sum(cart2theta,1);

figure;
imagesc(twoThetaAxis,chiAxis,cart2theta);
colormap('turbo')
caxis([0,150])
% daspect([1,1,1])

figure;
plot(twoThetaAxis,DiffProfile)
%% Integrate chi ranges
wedgeList = 1:Nwedges;
wedgeCenters = (window+1)/2:window:Nchi;
wedgeCentersDegree = wedgeCenters/(Nchi/360);
Wedges=struct;
intOffset = 0;
figure;
for ww = 1:Nwedges
    shiftedPattern = circshift(cart2theta,(1-ww)*window);
    Wedges(ww).wedgeCenter = wedgeCenters(ww);
    Wedges(ww).wedgeCenterDegree = wedgeCentersDegree(ww);
    Wedges(ww).DiffProfile = sum(shiftedPattern(1:window,:));
    plot(twoThetaAxis,Wedges(ww).DiffProfile+ww*intOffset)
    hold on
end
hold off

%%

% DiffProfile=sum(cart2theta,1);
% DiffProfile0deg=sum(cart2theta(1:floor((15/360)*Nchi),:),1);
% DiffProfile90deg=sum(cart2theta(floor((90/360)*Nchi):floor((105/360)*Nchi),:),1);
% DiffProfile180deg=sum(cart2theta(floor((180/360)*Nchi):floor((195/360)*Nchi),:),1);
% DiffProfile270deg=sum(cart2theta(floor((270/360)*Nchi):floor((285/360)*Nchi),:),1);
% 
% figure;
% plot(twoThetaAxis,DiffProfile0deg,twoThetaAxis,DiffProfile90deg,twoThetaAxis,DiffProfile180deg,twoThetaAxis,DiffProfile270deg,twoThetaAxis,DiffProfile*(15/360))
% legend('0 deg','90 deg','180 deg','270 deg','full')
% set(gca,'TickDir','out');

%% Fit around Quadruplet
xDiff=fit2thetaMaxInd-fit2thetaMinInd;
for ww = 1:Nwedges
    yDiff=Wedges(ww).DiffProfile(fit2thetaMaxInd)-Wedges(ww).DiffProfile(fit2thetaMinInd);
    newXaxis=[1:length(twoThetaAxis)]';
    newXaxis=newXaxis-fit2thetaMinInd;
    yAdjustment=(yDiff/xDiff)*newXaxis+Wedges(ww).DiffProfile(fit2thetaMinInd);
    AdjustedProfile=Wedges(ww).DiffProfile-yAdjustment';
    negativeIndices=find(AdjustedProfile<0); %find indices of all negative values
    AdjustedProfile(negativeIndices)=0; %set negative values to 0
    Wedges(ww).DiffProfileQuad=AdjustedProfile; % write to struct
end

for ww=1:Nwedges
    
    Wedges(ww).PeakFitIntensities2=zeros(num2Fit,1);
    Wedges(ww).PeakFitCenters2=zeros(num2Fit,1);
    Wedges(ww).PeakFitFWHMg2=zeros(num2Fit,1);
    Wedges(ww).PeakFitFWHMl2=zeros(num2Fit,1);
    Wedges(ww).PeakFitFWHM2=zeros(num2Fit,1);
    Wedges(ww).PeakFitEttas=zeros(num2Fit,1);    
        
    %perform fit for each wedge
    [tempFit,tempGOF]=fit(twoThetaAxis(fit2thetaMinInd:fit2thetaMaxInd)',Wedges(ww).DiffProfileQuad(fit2thetaMinInd:fit2thetaMaxInd)',PseudoVoigt4,'StartPoint',StartVals,'Upper',UpperBounds,'Lower',LowerBounds);

    %write fit information to persistent struct
    Wedges(ww).fitInfo=tempFit;
    Wedges(ww).GOF2=tempGOF;
    Wedges(ww).PeakFitIntensities2=[tempFit.a1,tempFit.a5,tempFit.a9,tempFit.a13];
    Wedges(ww).PeakFitCenters2=[tempFit.a2,tempFit.a6,tempFit.a10,tempFit.a14];
    Wedges(ww).PeakFitFWHMg2=[tempFit.a3,tempFit.a7,tempFit.a11,tempFit.a15];
    Wedges(ww).PeakFitFWHMl2=[tempFit.a4,tempFit.a8,tempFit.a12,tempFit.a16];
    Wedges(ww).PeakFitFWHM2=[FWHMcalc(tempFit.a3,tempFit.a4),FWHMcalc(tempFit.a7,tempFit.a8),FWHMcalc(tempFit.a11,tempFit.a12),FWHMcalc(tempFit.a15,tempFit.a16)];
    Wedges(ww).PeakFitEttas=[nCalc(tempFit.a3,tempFit.a4),nCalc(tempFit.a7,tempFit.a8),nCalc(tempFit.a11,tempFit.a12),nCalc(tempFit.a15,tempFit.a16)];    
    
end

%% Plot variation of peak centers with chi
kk = 1; %select peak to plot
PeakCenterList=zeros(1,Nwedges);
PeakIntensityList=zeros(1,Nwedges);
PeakFWHMList=zeros(1,Nwedges);
PeakGOFList=zeros(1,Nwedges);
for ww=1:Nwedges
    PeakCenterList(ww)=Wedges(ww).PeakFitCenters2(kk);
    PeakIntensityList(ww)=Wedges(ww).PeakFitIntensities2(kk);
    PeakFWHMList(ww)=Wedges(ww).PeakFitFWHM2(kk);
    PeakGOFList(ww)=Wedges(ww).GOF2.adjrsquare;
end

goodIndices = find(PeakGOFList>0.95);

figure;
plot(wedgeCentersDegree(goodIndices),PeakCenterList(goodIndices))
title(strcat('Peak-center of HAp ',PeakRef(kk),' vs. Azimuthal Angle'))
xlabel('Azimuthal Angle')
ylabel('Peak-center (2-theta)')
set(gca,'TickDir','out');

figure;
plot(wedgeCentersDegree(goodIndices),PeakIntensityList(goodIndices))
title(strcat('Peak Intensity of HAp ',PeakRef(kk),' vs. Azimuthal Angle'))
xlabel('Azimuthal Angle')
ylabel('Peak Intensity (a.u.)')
set(gca,'TickDir','out');

figure;
plot(wedgeCentersDegree(goodIndices),PeakFWHMList(goodIndices))
title(strcat('FWHM of HAp ',PeakRef(kk),' vs. Azimuthal Angle'))
xlabel('Azimuthal Angle')
ylabel('FWHM (2-theta)')
set(gca,'TickDir','out');

shiftedPeakCenterList = circshift(PeakCenterList,Nwedges/2);
figure;
plot(wedgeCentersDegree(goodIndices),PeakCenterList(goodIndices)-shiftedPeakCenterList(goodIndices))
title(strcat('Symmetric Difference in Peak Position of HAp ',PeakRef(kk),' vs. Azimuthal Angle'))
xlabel('Azimuthal Angle')
ylabel('Peak Center Difference (2-theta)')
set(gca,'TickDir','out');

figure;
plot(wedgeCentersDegree,PeakGOFList)
title(strcat('GOF of HAp ',PeakRef(kk),' vs. Azimuthal Angle'))
xlabel('Azimuthal Angle')
ylabel('Adjusted r-squared')
set(gca,'TickDir','out');

%% Fit azimuthal distribution with sinusoid to determine angle of beam center offset
Sinusoid1 = @(a1,a2,a3,x) a1*sind(x+a2)+a3;
StartVals2 = [0.05,50,0];
sinFit = fit(wedgeCentersDegree',(PeakCenterList-circshift(PeakCenterList,Nwedges/2))',Sinusoid1,'StartPoint',StartVals2)
figure;
plot(sinFit,wedgeCentersDegree,PeakCenterList-circshift(PeakCenterList,Nwedges/2))
xlabel('Azimuthal Angle')
ylabel('Peak Center Difference (2-theta)')
set(gca,'TickDir','out');

%% Determine refinements in x and y of beam center based on sinusoidal fit
pixelMag=radialPixelDistance(TwoThetaGuesses(kk)+abs(sinFit.a1)/2)-radialPixelDistance(TwoThetaGuesses(kk)-abs(sinFit.a1)/2); %convert the magnitude of the offset it 2-theta to pixels
shiftMagnitude=(sinFit.a1/abs(sinFit.a1))*pixelMag/2;
shiftPhase=-sinFit.a2+90;

xcNew=xc+shiftMagnitude*cosd(shiftPhase)
ycNew=yc+shiftMagnitude*sind(shiftPhase)

%% Main
%Loop over patterns in sample region
for j=pattSample
    %read in background corrected diffractin pattern .tif
    filename=strcat(strcat(scanName,'_',num2str(j,'%04.0f'),'.tif'));
    rawImage=double(imread(filename));    
    
    %create interpolant from image based on best guess for beam center
    F = griddedInterpolant(flipud(DY),DX,flipud(rawImage),'linear');
    
    %Convert from polar to cartesian coordinates (2theta)
    Pattern(j).Centers = struct;
    tempCenters = struct;
    for ii = 1:numCenters
        tempCenters(ii).xc=Xcenters(ii);
        tempCenters(ii).yc=Ycenters(ii);
        tempCenters(ii).queryPts = diffpol2cart_getQueryPts(twoThetaAxis,chiAxis,[Xoffsets(ii),Yoffsets(ii)],[tilt.rot,tilt.angle],Sdd);
        tempCenters(ii).cart2theta = reshape(F(tempCenters(ii).queryPts(2,:),tempCenters(ii).queryPts(1,:)),Nchi,NtwoTheta); 
    end
    
    %Compute average intensity (2theta) in each pattern
    centerIndex = ceil(numel(Xcenters)/2);
    Pattern(j).AvgInt2theta=sum(sum(tempCenters(centerIndex).cart2theta,1),2)/(NtwoTheta*Nchi);
    
    % Extract 1-D diffraction profiles from cartesian transforms: integrating 15 degree ranges at 0 90 180 and 270 degrees, as well as all 360 degrees
    for ii=1:numCenters
        tempCenters(ii).DiffProfile=sum(tempCenters(ii).cart2theta,1);
        tempCenters(ii).DiffProfile0deg=sum(tempCenters(ii).cart2theta(1:floor((15/360)*Nchi),:),1);
        tempCenters(ii).DiffProfile90deg=sum(tempCenters(ii).cart2theta(floor((90/360)*Nchi):floor((105/360)*Nchi),:),1);
        tempCenters(ii).DiffProfile180deg=sum(tempCenters(ii).cart2theta(floor((180/360)*Nchi):floor((195/360)*Nchi),:),1);
        tempCenters(ii).DiffProfile270deg=sum(tempCenters(ii).cart2theta(floor((270/360)*Nchi):floor((285/360)*Nchi),:),1);
    end    
    
    %Perform background subtraction around quadruplet
    for ii=1:numCenters
        xDiff=fit2thetaMaxInd-fit2thetaMinInd;
        yDiff=tempCenters(ii).DiffProfile(fit2thetaMaxInd)-tempCenters(ii).DiffProfile(fit2thetaMinInd);
        newXaxis=[1:length(twoThetaAxis)]';
        newXaxis=newXaxis-fit2thetaMinInd;
        yAdjustment=(yDiff/xDiff)*newXaxis+tempCenters(ii).DiffProfile(fit2thetaMinInd);
        AdjustedProfile=tempCenters(ii).DiffProfile-yAdjustment';
        negativeIndices=find(AdjustedProfile<0); %find indices of all negative values
        AdjustedProfile(negativeIndices)=0; %set negative values to 0
        tempCenters(ii).DiffProfileQuad=AdjustedProfile; % write to struct
    end
    
    %Perform fit of quadruplet for each center
    for ii=1:numCenters
        %initialize storage matrices
        Pattern(j).Centers(ii).PeakFitIntensities2=zeros(num2Fit,1);
        Pattern(j).Centers(ii).PeakFitCenters2=zeros(num2Fit,1);
        Pattern(j).Centers(ii).PeakFitFWHMg2=zeros(num2Fit,1);
        Pattern(j).Centers(ii).PeakFitFWHMl2=zeros(num2Fit,1);
        Pattern(j).Centers(ii).PeakFitFWHM2=zeros(num2Fit,1);
        Pattern(j).Centers(ii).PeakFitEttas=zeros(num2Fit,1);

        %perform fit for pattern in question
        [tempFit,tempGOF]=fit(twoThetaAxis(fit2thetaMinInd:fit2thetaMaxInd)',tempCenters(ii).DiffProfileQuad(fit2thetaMinInd:fit2thetaMaxInd)',PsuedoVoigt4,'StartPoint',StartVals,'Upper',UpperBounds,'Lower',LowerBounds);

        %write fit information to persistent struct
        Pattern(j).Centers(ii).fitInfo2=tempFit;
        Pattern(j).Centers(ii).GOF2=tempGOF;
        Pattern(j).Centers(ii).PeakFitIntensities2=[tempFit.a1,tempFit.a5,tempFit.a9,tempFit.a13];
        Pattern(j).Centers(ii).PeakFitCenters2=[tempFit.a2,tempFit.a6,tempFit.a10,tempFit.a14];
        Pattern(j).Centers(ii).PeakFitFWHMg2=[tempFit.a3,tempFit.a7,tempFit.a11,tempFit.a15];
        Pattern(j).Centers(ii).PeakFitFWHMl2=[tempFit.a4,tempFit.a8,tempFit.a12,tempFit.a16];
        Pattern(j).Centers(ii).PeakFitFWHM2=[FWHMcalc(tempFit.a3,tempFit.a4),FWHMcalc(tempFit.a7,tempFit.a8),FWHMcalc(tempFit.a11,tempFit.a12),FWHMcalc(tempFit.a15,tempFit.a16)];
        Pattern(j).Centers(ii).PeakFitEttas=[nCalc(tempFit.a3,tempFit.a4),nCalc(tempFit.a7,tempFit.a8),nCalc(tempFit.a11,tempFit.a12),nCalc(tempFit.a15,tempFit.a16)];
    end    
    
    %Compute crystallite size from FWHM information
    for ii=1:numCenters
        PeakFitBetas=zeros(num2Fit,1);
        PeakFitBetasCorrected=zeros(num2Fit,1);
        BetaStrain=(180/pi)*4*0.00055*tan((pi/360)*Pattern(j).Centers(ii).PeakFitCenters2); %Compute broadening due to average microstrain extracted from W-H of 0.00055. Bstrain = 4*epsilon*tan(theta) 

        %compute integral breadth (area under peak divided by max intensity)
        %for all reflections
        PeakFitBetas=Pattern(j).Centers(ii).PeakFitIntensities2./PsuedoVoigt(0,Pattern(j).Centers(ii).PeakFitIntensities2,0,Pattern(j).Centers(ii).PeakFitFWHMg2,Pattern(j).Centers(ii).PeakFitFWHMl2);
        smallInd=find(PeakFitBetas<SiFit.a1./PsuedoVoigt(0,SiFit.a1,0,SiFit.a3,SiFit.a4)); %determine if Beta is smaller than instrumental value from Si standard
        PeakFitBetasConstrained=PeakFitBetas;
        PeakFitBetasConstrained(smallInd)=(SiFit.a1./PsuedoVoigt(0,SiFit.a1,0,SiFit.a3,SiFit.a4))*1.000001; %make sure the value will still be positive after subtraction, but very small.

        %subtract instrumental broadening
        PeakFitBetasCorrectedSquared=PeakFitBetasConstrained.^2-(SiFit.a1./PsuedoVoigt(0,SiFit.a1,0,SiFit.a3,SiFit.a4)).^2;
        PeakFitBetasCorrected=sqrt(PeakFitBetasCorrectedSquared);

        %also subtract strain broadening computed through W-H analysis
        smallInd=find(PeakFitBetasCorrected<BetaStrain);
        PeakFitBetasCorrectedConstrained=PeakFitBetasCorrected;
        PeakFitBetasCorrectedConstrained(smallInd)=1.000001*BetaStrain(smallInd);
        PeakFitBetasCorrectedWH=PeakFitBetasCorrectedConstrained-BetaStrain;

        %Calculate Scherrer using K factors for specific reflections (being
        %careful to convert to radians) and save to persistent struct
        Pattern(j).Centers(ii).PeakScherrer=(lambda*Kfactors)./((pi/180)*PeakFitBetasCorrected.*cos((pi/360)*Pattern(j).Centers(ii).PeakFitCenters2));
        Pattern(j).Centers(ii).PeakScherrerWH=(lambda*Kfactors)./((pi/180)*PeakFitBetasCorrectedWH.*cos((pi/360)*Pattern(j).Centers(ii).PeakFitCenters2));
    end
    
    %Compute lattice parameter from mean peak position of 121 and 030
    %reflections
    for ii=1:numCenters
        d003 = lambda/(2*sind(Pattern(j).Centers(ii).PeakFitCenters2(3)/2));
        d121 = lambda/(2*sind(Pattern(j).Centers(ii).PeakFitCenters2(1)/2));
        Pattern(j).Centers(ii).a = 3*d003/cosd(30); %calculate a from 030 reflection
        Pattern(j).Centers(ii).c = ((1/d121^2)-((4/3)*7/(Pattern(j).Centers(ii).a^2)))^(-1/2); %calculate c from 121 and a parameter
        Pattern(j).Centers(ii).Vol=(sqrt(3)/2)*(Pattern(j).Centers(ii).a^2)*(Pattern(j).Centers(ii).c); %calculate volume for hexagonal crystal system
    end
    
    %Assemble Matrices to visualize across all patterns for each beam
    %center under analysis
    for ii=1:numCenters
        Centers(ii).AMatrix(ceil(j/xPts),(1+mod(j-1,xPts)))=Pattern(j).Centers(ii).a;
        Centers(ii).CMatrix(ceil(j/xPts),(1+mod(j-1,xPts)))=Pattern(j).Centers(ii).c;
        Centers(ii).VolMatrix(ceil(j/xPts),(1+mod(j-1,xPts)))=Pattern(j).Centers(ii).Vol;
        Centers(ii).Size121Matrix(ceil(j/xPts),(1+mod(j-1,xPts)))=Pattern(j).Centers(ii).PeakScherrerWH(1);
    end
    disp(j)
end
    
    
%% Visualize how the crystallographic parameters vary over beam center
CentersIndex = 1:numCenters;
CentersIndexMatrix = reshape(CentersIndex,[Ncenter_sampling,Ncenter_sampling])';

fig=figure(1);
set(fig,'Position',[200,200,1000,1000])
title('a parameter');
for ii=1:numCenters
    subplot(Ncenter_sampling,Ncenter_sampling,CentersIndexMatrix(ii))
    imagesc(Centers(ii).AMatrix(:,9:22))
    daspect([yAsp,xAsp,1])
    colormap(gca,parula);
%     title(strcat('a parameter - Center: ,', num2str(Xcenters(ii)),', ',num2str(Ycenters(ii))))
%     xlabel(XlabelStr)
%     ylabel(YlabelStr)
%     cbar = colorbar;
    caxis([9.438,9.453])
%     caxis('auto')
    set(gca,'TickDir','out');
%     cbar.TickDirection = 'out';
end

fig=figure(2);
set(fig,'Position',[200,200,1000,1000])
title('c parameter');
for ii=1:numCenters
    subplot(Ncenter_sampling,Ncenter_sampling,CentersIndexMatrix(ii))
    imagesc(Centers(ii).CMatrix(:,9:22))
    daspect([yAsp,xAsp,1])
    colormap(gca,parula);
%     title(strcat('c parameter - Center: ,', num2str(Xcenters(ii)),', ',num2str(Ycenters(ii))))
%     xlabel(XlabelStr)
%     ylabel(YlabelStr)
%     cbar = colorbar;
    caxis([6.85,6.91])
%     caxis('auto')
    set(gca,'TickDir','out');
%     cbar.TickDirection = 'out';
end

fig=figure(3);
set(fig,'Position',[200,200,1000,1000])
title('121 Size');
for ii=1:numCenters
    subplot(Ncenter_sampling,Ncenter_sampling,CentersIndexMatrix(ii))
    imagesc(Centers(ii).Size121Matrix(:,9:22))
    daspect([yAsp,xAsp,1])
    colormap(gca,parula);
%     title(strcat('Crystallite size || 121 - Center: ,', num2str(Xcenters(ii)),', ',num2str(Ycenters(ii))))
%     xlabel(XlabelStr)
%     ylabel(YlabelStr)
%     cbar = colorbar;
    caxis([230,330])
%     caxis('auto')
    set(gca,'TickDir','out');
%     cbar.TickDirection = 'out';
end

% fig=figure(4);
% set(fig,'Position',[200,200,1000,1000])
% for ii=1:numCenters
%     subplot(Ncenter_sampling,Ncenter_sampling,CentersIndexMatrix(ii))
%     imagesc(Centers(ii).VolMatrix(:,9:22))
%     daspect([yAsp,xAsp,1])
%     colormap(gca,parula);
%     title(strcat('Unit Cell Vol - Center: ,', num2str(Xcenters(ii)),', ',num2str(Ycenters(ii))))
%     xlabel(XlabelStr)
%     ylabel(YlabelStr)
%     cbar = colorbar;
%     caxis([530,533])
%     set(gca,'TickDir','out');
%     cbar.TickDirection = 'out';
% end
%% plot peak centers of 030 and 121 vs beam center position
peakCenterMatrix030 = zeros(size(Xcenters));
peakCenterMatrix121 = zeros(size(Xcenters));
for ii=1:numel(Xcenters)
    peakCenterMatrix121(ii)=tempCenters(ii).fitInfo2.a2;    
    peakCenterMatrix030(ii)=tempCenters(ii).fitInfo2.a10;
end

figure;
imagesc(xcenters_array,ycenters_array,peakCenterMatrix030)
title(strcat('Pattern ',num2str(jj),' 2 theta position of 030 relative to beam center'))
xlabel('beam center x')
ylabel('beam center y')
set(gca,'TickDir','out');

figure;
imagesc(xcenters_array,ycenters_array,peakCenterMatrix121)
title(strcat('Pattern ',num2str(jj),' 2 theta position of 121 relative to beam center'))
xlabel('beam center x')
ylabel('beam center y')
set(gca,'TickDir','out');


%% Plot a subset of fit datasets to check quality of fits
for j=pattSample(1:100:end)
    ii=5;
    xx=twoThetaAxis(fit2thetaMinInd:fit2thetaMaxInd);
    figure;
    plot(Pattern(j).Centers(ii).fitInfo2,'Blue')
    title(strcat('Pattern ',num2str(jj),', Center ',num2str(ii),' Fitting Check'))
    hold on
    plot(xx,Pattern(j).Centers(ii).DiffProfileQuad(fit2thetaMinInd:fit2thetaMaxInd),'o','Color','Black')
    plot(Pattern(j).Centers(ii).fitInfo2,xx',Pattern(j).Centers(ii).DiffProfileQuad(fit2thetaMinInd:fit2thetaMaxInd)','residuals')
    legend('fit','data','residual','zero line');
    hold off
    
    plot(Pattern(j).Centers(ii).fitInfo2,xx,Pattern(j).Centers(ii).DiffProfileQuad(fit2thetaMinInd:fit2thetaMaxInd))
    title(strcat('Pattern ',num2str(jj),', Center ',num2str(ii),' Fitting Check'))
    hold on
    plot(xx,PsuedoVoigt(xx,Pattern(j).Centers(ii).fitInfo2.a1,Pattern(j).Centers(ii).fitInfo2.a2,Pattern(j).Centers(ii).fitInfo2.a3,Pattern(j).Centers(ii).fitInfo2.a4),xx,PsuedoVoigt(xx,Pattern(j).Centers(ii).fitInfo2.a5,Pattern(j).Centers(ii).fitInfo2.a6,Pattern(j).Centers(ii).fitInfo2.a7,Pattern(j).Centers(ii).fitInfo2.a8),xx,PsuedoVoigt(xx,Pattern(j).Centers(ii).fitInfo2.a9,Pattern(j).Centers(ii).fitInfo2.a10,Pattern(j).Centers(ii).fitInfo2.a11,Pattern(j).Centers(ii).fitInfo2.a12),xx,PsuedoVoigt(xx,Pattern(j).Centers(ii).fitInfo2.a13,Pattern(j).Centers(ii).fitInfo2.a14,Pattern(j).Centers(ii).fitInfo2.a15,Pattern(j).Centers(ii).fitInfo2.a16))
    hold off
    
end





% clc%% Check beam calibration by comparing chi ranges from 0, 90, 180, and 270 degrees
% figure;
% jj=710;
% plot(twoThetaAxis,Pattern(jj).DiffProfile0deg,twoThetaAxis,Pattern(jj).DiffProfile90deg,twoThetaAxis,Pattern(jj).DiffProfile180deg,twoThetaAxis,Pattern(jj).DiffProfile270deg,twoThetaAxis,Pattern(jj).DiffProfile*(15/360))
% legend('0 deg','90 deg','180 deg','270 deg','full')
% set(gca,'TickDir','out');
% xlabel('1-theta')
% ylabel('Azimuthal angle from +x axis')



