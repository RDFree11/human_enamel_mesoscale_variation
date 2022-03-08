% Use this script with initial fit parameters generated in Fit2D to refine
% the beam center using the ceria standard. The code should be run iteratively
% to refine the beam center. 3 iterations are included in the current code.

% These analyses support the publication of:
% "Mesoscale Structure and Composition Varies Systematically in Human Tooth
% Enamel", by R. Free, K. DeRocher, V. Cooley, R. Xu, S.R. Stock, and D. Joester.
% This script also includes the units and axes information for each plot.

% Author: Robert Free
% Editor: Derk Joester

%% Clean Up
clear all;
close all;
clc
%% Flags and Paths
save_tables = true;
pn_out = './';
fn_out = 'TableS1.csv';
%% Initialize Structs to hold data
Centers=struct;
%% Define Parameters

%detector parameters from Fit2D and beamline
samp2det = 155.3818; %sample to detector distance in mm (may change per scan)
pixelSize = 0.079278; %pixel size in mm (usually the same for MAR detector from 34ide)
Sdd = samp2det / pixelSize; %sample to detector distance in pixels
% tilt.rot = 151.5569; % rotation of tilt plane relative to x-axis in degrees
% tilt.angle = 0.894298; % rotation of detector normal within the tilt plane in degrees

tilt.rot = 167.0; % rotation of tilt plane relative to x-axis in degrees
tilt.angle = 0.675; % rotation of detector normal within the tilt plane in degrees

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
energyKeV=17.000; %energy for scan
h=4.135667*10^-18; %keV*s
c=2.998792*10^18; %angstrom/s
lambda=h*c/energyKeV; %anstroms

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
filename='CeO2_17keV_dbl60s_001.tif';
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

%% Define beam center and beam center sampling matrices

% Beam center initializations for iterative refinement
% xc = 1050.303;
% yc = 887.4738;

%1st refinement
% xc = 1050.6435;
% yc = 888.39335;

%2nd refinement
% xc = 1050.614;
% yc = 888.315;

%zero tilt beam center for 111
% xc = 1049.073;
% yc = 889.1317;

%refinement update
% xc = xcNew;
% yc = ycNew;

%chosen center
xc = 1050.367;
yc = 888.856;

DX = X-xc;
DY = Y-yc;
%% Read in Ceria standard
rawImage=double(imread(filename)); 

figure;
imagesc(xAxis,yAxis,rawImage);
title('2D Ceria Diffraction Pattern')
caxis([0,200])
daspect([1,1,1])
xlabel('pixels')
ylabel('pixels')

F = griddedInterpolant(flipud(DY),DX,flipud(rawImage),'linear');
queryPts = diffpol2cart_getQueryPts(twoThetaAxis,chiAxis,[0,0],[tilt.rot,tilt.angle],Sdd);
cart2theta = reshape(F(queryPts(2,:),queryPts(1,:)),Nchi,NtwoTheta); 

figure;
imagesc(twoThetaAxis,chiAxis,cart2theta);
title('Cartesian Transform of Ceria Diffraction Pattern')
caxis([0,200])
% daspect([1,1,1])
xlabel('Bragg Angle (2-theta)')
ylabel('Azimuthal Angle (Chi)')
%% Integrate chi ranges
wedgeList = 1:Nwedges;
wedgeCenters = (window+1)/2:window:Nchi;
wedgeCentersDegree = wedgeCenters/(Nchi/360);
Wedges=struct;
intOffset = 500;
for ww = 1:Nwedges
    shiftedPattern = circshift(cart2theta,(1-ww)*window);
    Wedges(ww).wedgeCenter = wedgeCenters(ww);
    Wedges(ww).wedgeCenterDegree = wedgeCentersDegree(ww);
    Wedges(ww).DiffProfile = sum(shiftedPattern(1:window,:));
end

%% Define fit parameters (for 2theta fits of 111)

%specify initial values for fitting parameters (fit stability is sensitive
%to these. Determined through itereative fitting and checking quality of
%fit)
TwoThetaGuess=13.40; %two-theta guesses for reflection 220
IntGuess=8000; %intensity guess
FWHMgGuess=0.08; %FWHM 2theta of gaussian component guess
FWHMlGuess=0.08; %FWHM 2theta of lorentzian component guess
%combine all initial guesses into StartVals vector for fit function
StartVals=[IntGuess,TwoThetaGuess,FWHMgGuess,FWHMlGuess,]; %initial guesses for parameters

% Define limits for parameters to explore during fitting
twoThetaVariation=0.4; % max size of permissible variation in 2theta from initial guess
FWHMMin=min(SiFit.a3,SiFit.a4); %use the minimum from the silicon fit to determine minimum FWHM allowed during fitting.
FWHMMax=0.15; % max FWHM Variation (in either gaussian or lorentzian component)from initial guess allowed
IntMin=0.0001;
IntMax=10000000;
LowerBounds=[IntMin,TwoThetaGuess-twoThetaVariation,FWHMMin,FWHMMin];
UpperBounds=[IntMax,TwoThetaGuess+twoThetaVariation,FWHMMax,FWHMMax];
  %% Fit around 111 degree peak
fit2thetaMinInd=find(twoThetaAxis>13.2,1); %define index for left edge of target region (used for background subtraction)
fit2thetaMaxInd=find(twoThetaAxis<13.7,1,'last'); %and for right edge
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
    %perform fit for each wedge
    [tempFit,tempGOF]=fit(twoThetaAxis(fit2thetaMinInd:fit2thetaMaxInd)',Wedges(ww).DiffProfileQuad(fit2thetaMinInd:fit2thetaMaxInd)',PseudoVoigt1,'StartPoint',StartVals,'Upper',UpperBounds,'Lower',LowerBounds);

    %write fit information to persistent struct
    Wedges(ww).fitInfo=tempFit;
    Wedges(ww).GOF=tempGOF;
    Wedges(ww).PeakFitIntensity=tempFit.a1;
    Wedges(ww).PeakFitCenter=tempFit.a2;
    Wedges(ww).PeakFitFWHMg=tempFit.a3;
    Wedges(ww).PeakFitFWHMl=tempFit.a4;
    Wedges(ww).PeakFitFWHM=FWHMcalc(tempFit.a3,tempFit.a4);
    Wedges(ww).PeakFitEtta=nCalc(tempFit.a3,tempFit.a4);
end
    %% Plot variation of 111 peak center with chi
PeakCenterList=zeros(1,Nwedges);
PeakIntensityList=zeros(1,Nwedges);
PeakFWHMList=zeros(1,Nwedges);
for ww=1:Nwedges
    PeakCenterList(ww)=Wedges(ww).PeakFitCenter;
    PeakIntensityList(ww)=Wedges(ww).PeakFitIntensity;
    PeakFWHMList(ww)=Wedges(ww).PeakFitFWHM;
end

figure;
plot(wedgeCentersDegree,PeakCenterList)
title('Peak-center of Ceria 111 vs. Azimuthal Angle')
xlabel('Azimuthal Angle')
ylabel('Peak-center (2-theta)')
set(gca,'TickDir','out');

% figure;
% plot(wedgeCentersDegree,PeakIntensityList)
% title('Peak Intensity of Ceria 111 vs. Azimuthal Angle')
% xlabel('Azimuthal Angle')
% ylabel('Peak-center (2-theta)')
% set(gca,'TickDir','out');
% 
% figure;
% plot(wedgeCentersDegree,PeakFWHMList)
% title('FWHM of Ceria 111 vs. Azimuthal Angle')
% xlabel('Azimuthal Angle')
% ylabel('FWHM (2-theta)')
% set(gca,'TickDir','out');

figure;
plot(wedgeCentersDegree,PeakCenterList-mean(PeakCenterList))
title('Mean-difference of Peak Center of Ceria 111 vs. Azimuthal Angle')
xlabel('Azimuthal Angle')
ylabel('Peak-center (2-theta)')
set(gca,'TickDir','out');

figure;
plot(wedgeCentersDegree,PeakCenterList-circshift(PeakCenterList,Nwedges/2))
title('Symmetric Difference in Peak Position of Ceria 111 vs. Azimuthal Angle')
xlabel('Azimuthal Angle')
ylabel('Peak Center Difference (2-theta)')
set(gca,'TickDir','out');
  %% Fit azimuthal distribution with sinusoid to determine angle of beam center offset
Sinusoid1 = @(a1,a2,a3,x) a1*sind(x+a2)+a3;
StartVals2 = [0.05,50,0];
sinFit = fit((wedgeCentersDegree)',(PeakCenterList-circshift(PeakCenterList,Nwedges/2))',Sinusoid1,'StartPoint',StartVals2)
figure;
plot(sinFit,(wedgeCentersDegree),PeakCenterList-circshift(PeakCenterList,Nwedges/2))
xlabel('Azimuthal Angle')
ylabel('Peak Center Difference (2-theta)')
title('Symmetric Difference Fit With Sinusoid')
set(gca,'TickDir','out');
  %% Determine refinements in x and y of beam center based on sinusoidal fit (111)
pixelMag=radialPixelDistance(13.4+abs(sinFit.a1)/2)-radialPixelDistance(13.4-abs(sinFit.a1)/2); %convert the magnitude of the offset it 2-theta to pixels
shiftMagnitude=(sinFit.a1/abs(sinFit.a1))*pixelMag/2;
shiftPhase=-sinFit.a2+90;

xcNew=xc+shiftMagnitude*cosd(shiftPhase)
ycNew=yc+shiftMagnitude*sind(shiftPhase)
%% Define fit parameters (for 2theta fits of 220)

%specify initial values for fitting parameters (fit stability is sensitive
%to these. Determined through itereative fitting and checking quality of
%fit)
TwoThetaGuess=21.95; %two-theta guesses for reflection 220
IntGuess=8000; %intensity guess
FWHMgGuess=0.08; %FWHM 2theta of gaussian component guess
FWHMlGuess=0.08; %FWHM 2theta of lorentzian component guess
%combine all initial guesses into StartVals vector for fit function
StartVals=[IntGuess,TwoThetaGuess,FWHMgGuess,FWHMlGuess,]; %initial guesses for parameters

% Define limits for parameters to explore during fitting
twoThetaVariation=0.4; % max size of permissible variation in 2theta from initial guess
FWHMMin=min(SiFit.a3,SiFit.a4); %use the minimum from the silicon fit to determine minimum FWHM allowed during fitting.
FWHMMax=0.15; % max FWHM Variation (in either gaussian or lorentzian component)from initial guess allowed
IntMin=0.0001;
IntMax=10000000;
LowerBounds=[IntMin,TwoThetaGuess-twoThetaVariation,FWHMMin,FWHMMin];
UpperBounds=[IntMax,TwoThetaGuess+twoThetaVariation,FWHMMax,FWHMMax];
  %% Fit around 220 degree peak
fit2thetaMinInd=find(twoThetaAxis>21.4,1); %define index for left edge of target region (used for background subtraction)
fit2thetaMaxInd=find(twoThetaAxis<22.5,1,'last'); %and for right edge
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
    %perform fit for each wedge
    [tempFit,tempGOF]=fit(twoThetaAxis(fit2thetaMinInd:fit2thetaMaxInd)',Wedges(ww).DiffProfileQuad(fit2thetaMinInd:fit2thetaMaxInd)',PseudoVoigt1,'StartPoint',StartVals,'Upper',UpperBounds,'Lower',LowerBounds);

    %write fit information to persistent struct
    Wedges(ww).fitInfo=tempFit;
    Wedges(ww).GOF=tempGOF;
    Wedges(ww).PeakFitIntensity=tempFit.a1;
    Wedges(ww).PeakFitCenter=tempFit.a2;
    Wedges(ww).PeakFitFWHMg=tempFit.a3;
    Wedges(ww).PeakFitFWHMl=tempFit.a4;
    Wedges(ww).PeakFitFWHM=FWHMcalc(tempFit.a3,tempFit.a4);
    Wedges(ww).PeakFitEtta=nCalc(tempFit.a3,tempFit.a4);
end

  %% Plot variation of 221 peak center with chi
PeakCenterList=zeros(1,Nwedges);
PeakIntensityList=zeros(1,Nwedges);
PeakFWHMList=zeros(1,Nwedges);
for ww=1:Nwedges
    PeakCenterList(ww)=Wedges(ww).PeakFitCenter;
    PeakIntensityList(ww)=Wedges(ww).PeakFitIntensity;
    PeakFWHMList(ww)=Wedges(ww).PeakFitFWHM;
end

figure;
plot(wedgeCentersDegree,PeakCenterList)
title('Peak-center of Ceria 220 vs. Azimuthal Angle')
xlabel('Azimuthal Angle')
ylabel('Peak-center (2-theta)')
set(gca,'TickDir','out');

% figure;
% plot(wedgeCentersDegree,PeakIntensityList)
% title('Peak Intensity of Ceria 221 vs. Azimuthal Angle')
% xlabel('Azimuthal Angle')
% ylabel('Peak-center (2-theta)')
% set(gca,'TickDir','out');
% 
% figure;
% plot(wedgeCentersDegree,PeakFWHMList)
% title('FWHM of Ceria 221 vs. Azimuthal Angle')
% xlabel('Azimuthal Angle')
% ylabel('FWHM (2-theta)')
% set(gca,'TickDir','out');

figure;
plot(wedgeCentersDegree,PeakCenterList-mean(PeakCenterList))
title('Mean-difference of Peak Center of Ceria 220 vs. Azimuthal Angle')
xlabel('Azimuthal Angle')
ylabel('Peak-center (2-theta)')
set(gca,'TickDir','out');

figure;
plot(wedgeCentersDegree,PeakCenterList-circshift(PeakCenterList,Nwedges/2))
title('Symmetric Difference in Peak Position of Ceria 221 vs. Azimuthal Angle')
xlabel('Azimuthal Angle')
ylabel('Peak Center Difference (2-theta)')
set(gca,'TickDir','out');
  %% Fit azimuthal distribution with sinusoid to determine angle of beam center offset
Sinusoid1 = @(a1,a2,a3,x) a1*sind(x+a2)+a3;
StartVals2 = [0.05,50,0];
sinFit = fit((wedgeCentersDegree)',(PeakCenterList-circshift(PeakCenterList,Nwedges/2))',Sinusoid1,'StartPoint',StartVals2)
figure;
plot(sinFit,(wedgeCentersDegree),PeakCenterList-circshift(PeakCenterList,Nwedges/2))
xlabel('Azimuthal Angle')
ylabel('Peak Center Difference (2-theta)')
title('Symmetric Difference Fit With Sinusoid')
set(gca,'TickDir','out');
  %% Determine refinements in x and y of beam center based on sinusoidal fit (220)
pixelMag=radialPixelDistance(21.97+abs(sinFit.a1)/2)-radialPixelDistance(21.97-abs(sinFit.a1)/2); %convert the magnitude of the offset it 2-theta to pixels
shiftMagnitude=(sinFit.a1/abs(sinFit.a1))*pixelMag/2;
shiftPhase=-sinFit.a2+90;

xcNew=xc+shiftMagnitude*cosd(shiftPhase)
ycNew=yc+shiftMagnitude*sind(shiftPhase)
%% Output Table of Refined Parameters
Refinement={'Fit2D Calibrant','Matlab Refined'};
XrayEnergy={22.000,22.000};
XcenterPixel={1042.46,xcNew};
YcenterPixel={868.99,ycNew};
SampleToDetectorDistance={samp2det,samp2det};
TiltRotation={tilt.rot,tilt.rot};
TiltAngle={tilt.angle,tilt.angle};

sz = [2,7];
varTypes={'string','double','double','double','double','double','double'};
varNames={'Refinement','Xray Energy (kV)','X Center (pixel)','Y Center (pixel)','Sample-to-detector distance (mm)','Tilt Rotation (degrees)','Tilt Angle (degrees)'}

T1 = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
T1(:,1) = Refinement';
T1(:,2) = XrayEnergy';
T1(:,3) = XcenterPixel';
T1(:,4) = YcenterPixel';
T1(:,5) = SampleToDetectorDistance';
T1(:,6) = TiltRotation';
T1(:,7) = TiltAngle';

if save_tables
    writetable(T1,[pn_out,fn_out]);    
else
    disp('Table not saved!');
end