% Use this script to examine the robustness of the analysis to small errors
% in the assignment of the beam center. Currently, deviations of +/- 0.1
% pixels in x and y are examined.

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

%% Flags
plot_figures = false;
save_maps = false;
%% Initialize Strcuts to hold data
Pattern=struct;
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

%Define dimensions of map in pixels
xPts=23;
yPts=41;

%number of patterns
numPatt=xPts*yPts;
pattRange=1:numPatt;
pattMatrix=reshape(pattRange,[xPts,yPts])';
pattSampleMatrix=pattMatrix(:,9:22); %Cropped area based on examining results
pattSample=reshape(pattSampleMatrix',[1,numel(pattSampleMatrix)]); %Create list of patterns determined to be sample
num2Fit=4; %reflections in quadruplet

%experimental parameters
energyKeV=17.00000; %energy for scan
h=4.135667*10^-18; %keV*s
c=2.998792*10^18; %angstrom/s
lambda=h*c/energyKeV; %anstroms

%detector parameters
samp2det = 109.235; %sample to detector distance in mm
pixelSize = 0.079278; %pixel size in mm
Sdd = samp2det / pixelSize; %sample to detector distance in pixels
tilt.rot = 14.662552; % rotation of tilt plane relative to x-axis in degrees
tilt.angle = -1.327451; % rotation of detector normal within the tilt plane in degrees

%resolution definitions for 2-theta and Chi spaces
NtwoTheta = 1482;
Nchi = 720;
twoThetaAxis = linspace(atand(inner_lim/Sdd),atand(outer_lim/Sdd),NtwoTheta);
chiAxis = linspace(0,360,Nchi+1); 
chiAxis = chiAxis(1:end-1);
Ntwotheta = numel(twoThetaAxis);

Kfactors=[1,1,1,1]; %Define Scherrer factors for each reflection

%functions relating angular width on detector at different radial positions
delta2theta = @(radDis) atand((radDis+(pixelSize/2))/samp2det)-atand((radDis-(pixelSize/2))/samp2det);
radialPixelDistance = @(twoTheta) Sdd*tand(twoTheta);

% define name of sample and scan (important for reading data)
sampleName='Enamel D, ';
scanName='50s_pattern';
%% Define fitting functions (check working directory for referenced functions)

%set up combination of 4 Pseudo-voigt profiles as anonymous function
PseudoVoigt4 = @(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,x) PseudoVoigt(x,a1,a2,a3,a4)+PseudoVoigt(x,a5,a6,a7,a8)+PseudoVoigt(x,a9,a10,a11,a12)+PseudoVoigt(x,a13,a14,a15,a16);
PseudoVoigt2 = @(a1,a2,a3,a4,a5,a6,a7,a8,x) PseudoVoigt(x,a1,a2,a3,a4)+PseudoVoigt(x,a5,a6,a7,a8);
PseudoVoigt1 = @(a1,a2,a3,a4,x) PseudoVoigt(x,a1,a2,a3,a4);

%define combined full width half max function for Pseudo-voigt profile
FWHMcalc = @(fg,fl) ((fg.^5)+2.69269.*(fg.^4).*(fl)+2.42843.*(fg.^3).*(fl.^2)+4.47163.*(fg.^2).*(fl.^3)+0.07842.*(fg).*(fl.^4)+(fl.^5)).^(1/5);

%define etta computation for Pseudo voigt
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
twoThetaVariation=0.075; % max size of permissible variation in 2theta from initial guess
FWHMMin=min(SiFit.a3,SiFit.a4); %use the minimum from the silicon fit to determine minimum FWHM allowed during fitting.
FWHMMax=0.15; % max FWHM Variation (in either gaussian or lorentzian component)from initial guess allowed
IntMin=0.0001;
IntMax=10000000;
LowerBounds=[IntMin,TwoThetaGuesses(1)-twoThetaVariation,FWHMMin,FWHMMin,IntMin,TwoThetaGuesses(2)-twoThetaVariation,FWHMMin,FWHMMin,IntMin,TwoThetaGuesses(3)-twoThetaVariation,FWHMMin,FWHMMin,IntMin,TwoThetaGuesses(4)-twoThetaVariation,FWHMMin,FWHMMin];
UpperBounds=[IntMax,TwoThetaGuesses(1)+twoThetaVariation,FWHMMax,FWHMMax,IntMax,TwoThetaGuesses(2)+twoThetaVariation,FWHMMax,FWHMMax,IntMax,TwoThetaGuesses(3)+twoThetaVariation,FWHMMax,FWHMMax,IntMax,TwoThetaGuesses(4)+twoThetaVariation,FWHMMax,FWHMMax];

fit2thetaMinInd=find(twoThetaAxis>14.0,1); %define index for left edge of quadruplet region (used for background subtraction)
fit2thetaMaxInd=find(twoThetaAxis<16.3,1,'last'); %and for right edge
%% Define beam center and beam center sampling matrices

xc = 1041.90; %refined fit of beam center in x direction from ceria (indexed from left side of pattern)
yc = 869.3162; %refined fit of beam center in y direction from ceria (indexed from top of pattern)

Ncenter_sampling = 3; %number of samples (should be odd to include current beam center)
center_range = 0.1; %define size of steps of sampling in pixels

xoffset_array = linspace(-center_range,center_range,Ncenter_sampling);
yoffset_array = linspace(-center_range,center_range,Ncenter_sampling);

xcenters_array = xoffset_array+xc;
ycenters_array = yoffset_array+yc;

[Xoffsets,Yoffsets]=meshgrid(xoffset_array,yoffset_array);
[Xcenters,Ycenters]=meshgrid(xcenters_array,ycenters_array);

DX = X-xc;
DY = Y-yc;

numCenters = numel(Xcenters);
%initialize persistent struct to hold full sample maps of key parameters
for ii=1:numCenters
    Centers(ii).AMatrix=zeros(yPts,xPts);
    Centers(ii).CMatrix=zeros(yPts,xPts);
    Centers(ii).VolMatrix=zeros(yPts,xPts);
    Centers(ii).Size121Matrix=zeros(yPts,xPts);
    Centers(ii).FWHMMatrix=zeros(yPts,xPts);
end

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
        [tempFit,tempGOF]=fit(twoThetaAxis(fit2thetaMinInd:fit2thetaMaxInd)',tempCenters(ii).DiffProfileQuad(fit2thetaMinInd:fit2thetaMaxInd)',PseudoVoigt4,'StartPoint',StartVals,'Upper',UpperBounds,'Lower',LowerBounds);

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
        PeakFitBetas=Pattern(j).Centers(ii).PeakFitIntensities2./PseudoVoigt(0,Pattern(j).Centers(ii).PeakFitIntensities2,0,Pattern(j).Centers(ii).PeakFitFWHMg2,Pattern(j).Centers(ii).PeakFitFWHMl2);
        smallInd=find(PeakFitBetas<SiFit.a1./PseudoVoigt(0,SiFit.a1,0,SiFit.a3,SiFit.a4)); %determine if Beta is smaller than instrumental value from Si standard
        PeakFitBetasConstrained=PeakFitBetas;
        PeakFitBetasConstrained(smallInd)=(SiFit.a1./PseudoVoigt(0,SiFit.a1,0,SiFit.a3,SiFit.a4))*1.000001; %make sure the value will still be positive after subtraction, but very small.

        %subtract instrumental broadening
        PeakFitBetasCorrectedSquared=PeakFitBetasConstrained.^2-(SiFit.a1./PseudoVoigt(0,SiFit.a1,0,SiFit.a3,SiFit.a4)).^2;
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

if plot_figures
    fig=figure('Name','a parameter');
    set(fig,'Position',[200,200,1000,1000])
    for ii=1:numCenters
        subplot(Ncenter_sampling,Ncenter_sampling,CentersIndexMatrix(ii))
        imagesc(Centers(ii).AMatrix(:,9:22))
        daspect([yAsp,xAsp,1])
        colormap(gca,parula);
        title(strcat(num2str(Xcenters(ii)),', ',num2str(Ycenters(ii))))
        caxis([9.438,9.448])
        set(gca,'TickDir','out');
    end

    fig=figure('Name','c parameter');
    set(fig,'Position',[200,200,1000,1000])
    for ii=1:numCenters
        subplot(Ncenter_sampling,Ncenter_sampling,CentersIndexMatrix(ii))
        imagesc(Centers(ii).CMatrix(:,9:22))
        daspect([yAsp,xAsp,1])
        colormap(gca,parula);
        title(strcat(num2str(Xcenters(ii)),', ',num2str(Ycenters(ii))))
        caxis([6.85,6.92])
        set(gca,'TickDir','out');
    end

    fig=figure('Name','121 Size');
    set(fig,'Position',[200,200,1000,1000])
    for ii=1:numCenters
        subplot(Ncenter_sampling,Ncenter_sampling,CentersIndexMatrix(ii))
        imagesc(Centers(ii).Size121Matrix(:,9:22))
        daspect([yAsp,xAsp,1])
        colormap(gca,parula);
        title(strcat(num2str(Xcenters(ii)),', ',num2str(Ycenters(ii))))
        caxis([230,330])
        set(gca,'TickDir','out');
    end
end

%% Save Maps

if save_maps
    for ii=1:numCenters
        writematrix(Centers(ii).Size121Matrix(:,9:22),['121size_BCvariance_',num2str(Centers(ii))]);
        writematrix(Centers(ii).AMatrix(:,9:22),['aParameter_BCvariance_',num2str(Centers(ii))]);
        writematrix(Centers(ii).CMatrix(:,9:22),['cParameter_BCvariance_',num2str(Centers(ii))]);
    end
else
    disp('Maps Not Saved!')
end