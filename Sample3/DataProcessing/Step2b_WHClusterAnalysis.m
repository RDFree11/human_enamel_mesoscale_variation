% Use this script with the output of Step1_ProcessEnamelDiffractionPatterns
% and Step2_ClusterAnalysis to then conduct Williamson Hall analyses of
% composite whole sample and clustered diffraction patterns.

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
plot_patterns = true;
save_patterns = false;
save_results = false;
plot_WH = true;
%% Initialize Structs to hold data
Cluster=struct;
%% Define Parameters

%define detector axes
xAxis = linspace(1,2048,2048);
yAxis = linspace(2048,1,2048);
[Y,X] = ndgrid(yAxis,xAxis);
inner_lim = 50; %inner pixel limit for transform plot relative to beam center
outer_lim = 850; %outer pixel limit for transform plot relative to beam center

%define aspect ratio for plots and commonly used items for plots
xAsp=1;
yAsp=1;
XlabelStr='X-position (0.5um steps)';
YlabelStr='Y-position (0.5um steps)';
figWidth=300;
figHeight=500;

%Define dimensions of map in pixels (different for each scan)
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
samp2det = 109.2134; %sample to detector distance in mm (may change per scan)
pixelSize = 0.079278; %pixel size in mm (usually the same for MAR detector from 34ide)
Sdd = samp2det / pixelSize; %sample to detector distance in pixels
tilt.rot = 193.3; % rotation of tilt plane relative to x-axis in degrees
tilt.angle = 1.33; % rotation of detector normal within the tilt plane in degrees
% xc = 1041.940; %refined fit of beam center in x direction from ceria (indexed from left side of pattern)
% yc = 869.356; %refined fit of beam center in y direction from ceria (indexed from top of pattern)
xc = 1041.924; %refined fit from fit of 121 reflection of pattern 710
yc = 869.520; %refined fit from fit of 121 reflection of pattern 710

DX = X-xc;
DY = Y-yc;

NtwoTheta = 1482;
Nchi = 720;
rRes = NtwoTheta;
chiRes = NtwoTheta;
twoThetaAxis = linspace(atand(inner_lim/Sdd),atand(outer_lim/Sdd),NtwoTheta);
chiAxis = linspace(0,360,Nchi+1); 
chiAxis = chiAxis(1:end-1);
Ntwotheta = numel(twoThetaAxis);

Kfactors=[1,1,1,1]; %Define Scherrer factors for each reflection

%functions relating angular width on detector at different radial positions
delta2theta = @(radDis) atand((radDis+(pixelSize/2))/samp2det)-atand((radDis-(pixelSize/2))/samp2det);
radialPixelDistance = @(twoTheta) Sdd*tand(twoTheta);

% define name of sample and scan (important for reading data)
folderName = 'TotalMapCorrected/';
sampleName='Sample 3';
scanName='50s_pattern';

%% Read in cluster map assignments
% import cluster map and identify nan patterns
G  = readmatrix('Results 2 cluster/Cluster_Map.csv');
Gp = G(:);
nanPatt = find(isnan(Gp));
GpClean=Gp;
GpClean(nanPatt,:)=[];

%% Integrate Patterns together by cluster
RHIntMatrix=zeros(2048);
IRIntMatrix=zeros(2048);
SampleIntMatrix=zeros(2048);
%Loop over patterns in definied sample region to integrate patterns into
%either RH or IR/RT combined pattern based on cluster assignment
for j=1:numel(G)
    pattInd=pattSampleMatrix(j);
    clusterID=G(j);
    
    %read in background corrected diffractin pattern .tif
    filename=strcat(strcat(folderName,scanName,'_',num2str(pattInd,'%04.0f'),'.tif'));
    rawImage=double(imread(filename));

    SampleIntMatrix=SampleIntMatrix+rawImage;
    if clusterID==1
        RHIntMatrix=RHIntMatrix+rawImage;
    elseif clusterID==2
        IRIntMatrix=IRIntMatrix+rawImage;
    end
    disp(j)
end
    
%% Convert integrated patterns to cartesian coordinates

Cluster(1).name = 'Rod Head';
Cluster(2).name = 'Interrod/Rod Tail';
Cluster(3).name = 'Whole Sample';
%create interpolant from image based on best guess for beam center
F_RH = griddedInterpolant(flipud(DY),DX,flipud(RHIntMatrix),'linear');
F_IR = griddedInterpolant(flipud(DY),DX,flipud(IRIntMatrix),'linear');
F_samp = griddedInterpolant(flipud(DY),DX,flipud(SampleIntMatrix),'linear');

%Convert from polar to cartesian coordinates (2theta)
queryPts = diffpol2cart_getQueryPts(twoThetaAxis,chiAxis,[0,0],[tilt.rot,tilt.angle],Sdd);
Cluster(1).cart2theta = reshape(F_RH(queryPts(2,:),queryPts(1,:)),Nchi,NtwoTheta);
Cluster(2).cart2theta = reshape(F_IR(queryPts(2,:),queryPts(1,:)),Nchi,NtwoTheta);
Cluster(3).cart2theta = reshape(F_samp(queryPts(2,:),queryPts(1,:)),Nchi,NtwoTheta);

%Compute average intensity (2theta) in each pattern
Cluster(1).AvgInt2theta=sum(sum(Cluster(1).cart2theta,1),2)/(NtwoTheta*Nchi);
Cluster(2).AvgInt2theta=sum(sum(Cluster(2).cart2theta,1),2)/(NtwoTheta*Nchi);
Cluster(3).AvgInt2theta=sum(sum(Cluster(3).cart2theta,1),2)/(NtwoTheta*Nchi);

%Compute 1D diffraction profile vs 2-theta
Cluster(1).DiffProfile=sum(Cluster(1).cart2theta,1);
Cluster(2).DiffProfile=sum(Cluster(2).cart2theta,1);
Cluster(3).DiffProfile=sum(Cluster(3).cart2theta,1);

if plot_patterns
    close all
    figure
    plot(twoThetaAxis,Cluster(1).DiffProfile)
    figure
    plot(twoThetaAxis,Cluster(2).DiffProfile)
    figure
    plot(twoThetaAxis,Cluster(3).DiffProfile)
end

%% Define fitting functions (check working directory for referenced functions)

%set up combination of 4 pseudo-voigt profiles as anonymous function
PseudoVoigt4 = @(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,x) PseudoVoigt(x,a1,a2,a3,a4)+PseudoVoigt(x,a5,a6,a7,a8)+PseudoVoigt(x,a9,a10,a11,a12)+PseudoVoigt(x,a13,a14,a15,a16);
PseudoVoigt2 = @(a1,a2,a3,a4,a5,a6,a7,a8,x) PseudoVoigt(x,a1,a2,a3,a4)+PseudoVoigt(x,a5,a6,a7,a8);
PseudoVoigt1 = @(a1,a2,a3,a4,x) PseudoVoigt(x,a1,a2,a3,a4);

%define combined full width half max function for pseudo-voigt profile
FWHMcalc = @(fg,fl) ((fg.^5)+2.69269.*(fg.^4).*(fl)+2.42843.*(fg.^3).*(fl.^2)+4.47163.*(fg.^2).*(fl.^3)+0.07842.*(fg).*(fl.^4)+(fl.^5)).^(1/5);

%define etta computation for pseudo voigt
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
%% Williamson-Hall Analysis and Plotting
% close all

%define peak centers to target for Williamson-Hall plot
peakLabels={'010','020','1-30','1-31','030','031','130','222','2-52','240 and 331'};
peakLabels2={'010','020','1-30','1-31','030','031','130','222','2-52','240', '331'};
peakCenters=[5.0978,10.2344,13.5453,14.8603,15.3709,16.5467,18.5116,21.5904,25.5821,27.3947]; %centers of ranges to fit over (effectively peak centers except for last entry, which is average of two reflections)
FitRange=20; %number of indices to span to the right and left of peak for fitting

%specify initial values for fitting parameters
TwoThetaGuesses=[5.0978,10.2344,13.5453,14.8603,15.3709,16.5467,18.5116,21.5904,25.5821,27.2994,27.4886]; %peak centers of all actual peaks
IntGuess=500000; %intensity guess
FWHMgGuess=0.08; %FWHM 2theta of gaussian component guess
FWHMlGuess=0.08; %FWHM 2theta of lorentzian component guess
%combine all initial guesses into StartVals vector for fit function

% Define limits for parameters to explore during fitting
twoThetaVariation=0.075; % max size of permissible variation in 2theta from initial guess
FWHMMin=min(SiFit.a3,SiFit.a4); %use the minimum from the silicon fit to determine minimum FWHM allowed during fitting.
FWHMMax=0.15; % max FWHM Variation (in either gaussian or lorentzian component)from initial guess allowed
IntMin=0.0001;
IntMax=100000000;

for j=1:3    
    Cluster(j).PeakBetas=zeros(length(TwoThetaGuesses),1);
    Cluster(j).PeakBetasCorrected=zeros(length(TwoThetaGuesses),1);
    Cluster(j).PeakCenters=zeros(length(TwoThetaGuesses),1);
    Cluster(j).PeakFWHMs=zeros(length(TwoThetaGuesses),1);
    Cluster(j).PeakFWHMsCorrected=zeros(length(TwoThetaGuesses),1);
    Cluster(j).PeakFWHMs2=zeros(length(TwoThetaGuesses),1);
    Cluster(j).PeakFWHMs2Corrected=zeros(length(TwoThetaGuesses),1);

    for k=1:length(peakCenters)
   
        if k<10 %perform single background subtraction and fits
            centerIndex=find(twoThetaAxis>peakCenters(k),1);
            WHRangeIndices=(centerIndex-FitRange):(centerIndex+FitRange);
    
            Cluster(j).peak(k).label=peakLabels(k);
            Cluster(j).peak(k).twoThetaDomain=twoThetaAxis(WHRangeIndices);
            Cluster(j).peak(k).IntensityProfile=Cluster(j).DiffProfile(WHRangeIndices);
    
            %background subtraction for peak k of pattern j
            fit2thetaMinInd=centerIndex-FitRange; %define index for left edge of region to background subtract
            fit2thetaMaxInd=centerIndex+FitRange; %and for right edge
            xDiff=fit2thetaMaxInd-fit2thetaMinInd;
            yDiff=Cluster(j).DiffProfile(fit2thetaMaxInd)-Cluster(j).DiffProfile(fit2thetaMinInd);
            newXaxis=[1:length(twoThetaAxis)]';
            newXaxis=newXaxis-fit2thetaMinInd;
            yAdjustment=(yDiff/xDiff)*newXaxis+Cluster(j).DiffProfile(fit2thetaMinInd);
            AdjustedProfile=Cluster(j).DiffProfile-yAdjustment';
            negativeIndices=find(AdjustedProfile<0); %find indices of all negative values
            AdjustedProfile(negativeIndices)=0; %set negative values to 0
            Cluster(j).peak(k).AdjIntensityProfile=AdjustedProfile; % write to struct for peak k
    
            StartVals=[IntGuess,TwoThetaGuesses(k),FWHMgGuess,FWHMlGuess];
            LowerBounds=[IntMin,TwoThetaGuesses(k)-twoThetaVariation,FWHMMin,FWHMMin];
            UpperBounds=[IntMax,TwoThetaGuesses(k)+twoThetaVariation,FWHMMax,FWHMMax];
    
            [tempFit,tempGOF]=fit(twoThetaAxis(fit2thetaMinInd:fit2thetaMaxInd)',Cluster(j).peak(k).AdjIntensityProfile(fit2thetaMinInd:fit2thetaMaxInd)',PseudoVoigt1,'StartPoint',StartVals,'Upper',UpperBounds,'Lower',LowerBounds);
    
            xxx = -2:0.001:2;
            FitMax=PseudoVoigt(0,tempFit.a1,0,tempFit.a3,tempFit.a4);
            leftEdgeIndex = find(PseudoVoigt(xxx,tempFit.a1,0,tempFit.a3,tempFit.a4)>(FitMax/2),1);
            rightEdgeIndex = find(PseudoVoigt(xxx,tempFit.a1,0,tempFit.a3,tempFit.a4)>(FitMax/2),1,'last');
            Cluster(j).peak(k).PeakFitFWHM2 = xxx(rightEdgeIndex)-xxx(leftEdgeIndex);
            Cluster(j).PeakFWHMs2(k) = xxx(rightEdgeIndex)-xxx(leftEdgeIndex);
    
            Cluster(j).peak(k).fitInfo=tempFit;
            Cluster(j).peak(k).GOF=tempGOF;
            Cluster(j).peak(k).PeakFitIntensity=tempFit.a1;
            Cluster(j).peak(k).PeakFitCenter=tempFit.a2;
            Cluster(j).peak(k).PeakFitFWHMg=tempFit.a3;
            Cluster(j).peak(k).PeakFitFWHMl=tempFit.a4;
            Cluster(j).peak(k).PeakFitFWHM=FWHMcalc(tempFit.a3,tempFit.a4);
            Cluster(j).PeakFWHMs(k)=FWHMcalc(tempFit.a3,tempFit.a4);
    
            Cluster(j).peak(k).Beta=Cluster(j).peak(k).PeakFitIntensity./PseudoVoigt(0,Cluster(j).peak(k).PeakFitIntensity,0,Cluster(j).peak(k).PeakFitFWHMg,Cluster(j).peak(k).PeakFitFWHMl);
            Cluster(j).PeakBetas(k)=Cluster(j).peak(k).Beta;
            Cluster(j).PeakCenters(k)=Cluster(j).peak(k).PeakFitCenter;
        
        else %perform background subtraction and double fit for the 240 and 331

            centerIndex=find(twoThetaAxis>peakCenters(k),1);
            WHRangeIndices=(centerIndex-FitRange):(centerIndex+FitRange);
    
            Cluster(j).peak(k).label=peakLabels(k);
            Cluster(j).peak(k).twoThetaDomain=twoThetaAxis(WHRangeIndices);
            Cluster(j).peak(k).IntensityProfile=Cluster(j).DiffProfile(WHRangeIndices);
    
            %background subtraction for peak k of pattern j
            fit2thetaMinInd=centerIndex-FitRange; %define index for left edge of region to background subtract
            fit2thetaMaxInd=centerIndex+FitRange; %and for right edge
            xDiff=fit2thetaMaxInd-fit2thetaMinInd;
            yDiff=Cluster(j).DiffProfile(fit2thetaMaxInd)-Cluster(j).DiffProfile(fit2thetaMinInd);
            newXaxis=[1:length(twoThetaAxis)]';
            newXaxis=newXaxis-fit2thetaMinInd;
            yAdjustment=(yDiff/xDiff)*newXaxis+Cluster(j).DiffProfile(fit2thetaMinInd);
            AdjustedProfile=Cluster(j).DiffProfile-yAdjustment';
            negativeIndices=find(AdjustedProfile<0); %find indices of all negative values
            AdjustedProfile(negativeIndices)=0; %set negative values to 0
            Cluster(j).peak(k).AdjIntensityProfile=AdjustedProfile; % write to struct for peak k
    
            StartVals=[IntGuess,TwoThetaGuesses(k),FWHMgGuess,FWHMlGuess,IntGuess,TwoThetaGuesses(k+1),FWHMgGuess,FWHMlGuess];
            LowerBounds=[IntMin,TwoThetaGuesses(k)-twoThetaVariation,FWHMMin,FWHMMin,IntMin,TwoThetaGuesses(k+1)-twoThetaVariation,FWHMMin,FWHMMin];
            UpperBounds=[IntMax,TwoThetaGuesses(k)+twoThetaVariation,FWHMMax,FWHMMax,IntMax,TwoThetaGuesses(k+1)+twoThetaVariation,FWHMMax,FWHMMax];
    
            [tempFit,tempGOF]=fit(twoThetaAxis(fit2thetaMinInd:fit2thetaMaxInd)',Cluster(j).peak(k).AdjIntensityProfile(fit2thetaMinInd:fit2thetaMaxInd)',PseudoVoigt2,'StartPoint',StartVals,'Upper',UpperBounds,'Lower',LowerBounds);
    
            xxx = -2:0.001:2;
            FitMax=PseudoVoigt(0,tempFit.a1,0,tempFit.a3,tempFit.a4);
            leftEdgeIndex = find(PseudoVoigt(xxx,tempFit.a1,0,tempFit.a3,tempFit.a4)>(FitMax/2),1);
            rightEdgeIndex = find(PseudoVoigt(xxx,tempFit.a1,0,tempFit.a3,tempFit.a4)>(FitMax/2),1,'last');
            Cluster(j).peak(k).PeakFitFWHM2 = xxx(rightEdgeIndex)-xxx(leftEdgeIndex);
            Cluster(j).PeakFWHMs2(k) = xxx(rightEdgeIndex)-xxx(leftEdgeIndex);
    
            xxx = -2:0.001:2;
            FitMax=PseudoVoigt(0,tempFit.a5,0,tempFit.a7,tempFit.a8);
            leftEdgeIndex = find(PseudoVoigt(xxx,tempFit.a5,0,tempFit.a7,tempFit.a8)>(FitMax/2),1);
            rightEdgeIndex = find(PseudoVoigt(xxx,tempFit.a5,0,tempFit.a7,tempFit.a8)>(FitMax/2),1,'last');
            Cluster(j).peak(k+1).PeakFitFWHM2 = xxx(rightEdgeIndex)-xxx(leftEdgeIndex);
            Cluster(j).PeakFWHMs2(k+1) = xxx(rightEdgeIndex)-xxx(leftEdgeIndex);
    
            Cluster(j).peak(k).fitInfo=tempFit;
            Cluster(j).peak(k).GOF=tempGOF;
            Cluster(j).peak(k).PeakFitIntensity=tempFit.a1;
            Cluster(j).peak(k).PeakFitCenter=tempFit.a2;
            Cluster(j).peak(k).PeakFitFWHMg=tempFit.a3;
            Cluster(j).peak(k).PeakFitFWHMl=tempFit.a4;
            Cluster(j).peak(k).PeakFitFWHM=FWHMcalc(tempFit.a3,tempFit.a4);
            Cluster(j).PeakFWHMs(k)=FWHMcalc(tempFit.a3,tempFit.a4);
    
            Cluster(j).peak(k+1).fitInfo=tempFit;
            Cluster(j).peak(k+1).GOF=tempGOF;
            Cluster(j).peak(k+1).PeakFitIntensity=tempFit.a5;
            Cluster(j).peak(k+1).PeakFitCenter=tempFit.a6;
            Cluster(j).peak(k+1).PeakFitFWHMg=tempFit.a7;
            Cluster(j).peak(k+1).PeakFitFWHMl=tempFit.a8;
            Cluster(j).peak(k+1).PeakFitFWHM=FWHMcalc(tempFit.a7,tempFit.a8);
            Cluster(j).PeakFWHMs(k+1)=FWHMcalc(tempFit.a7,tempFit.a8);
    
            Cluster(j).peak(k).Beta=Cluster(j).peak(k).PeakFitIntensity./PseudoVoigt(0,Cluster(j).peak(k).PeakFitIntensity,0,Cluster(j).peak(k).PeakFitFWHMg,Cluster(j).peak(k).PeakFitFWHMl);
            Cluster(j).PeakBetas(k)=Cluster(j).peak(k).Beta;
            Cluster(j).peak(k+1).Beta=Cluster(j).peak(k+1).PeakFitIntensity./PseudoVoigt(0,Cluster(j).peak(k+1).PeakFitIntensity,0,Cluster(j).peak(k+1).PeakFitFWHMg,Cluster(j).peak(k+1).PeakFitFWHMl);
            Cluster(j).PeakBetas(k+1)=Cluster(j).peak(k+1).Beta;
            Cluster(j).PeakCenters(k)=Cluster(j).peak(k).PeakFitCenter;
            Cluster(j).PeakCenters(k+1)=Cluster(j).peak(k+1).PeakFitCenter;
        end
    end

    Cluster(j).PeakBetasCorrectedSquared=Cluster(j).PeakBetas.^2-(SiFit.a1./PseudoVoigt(0,SiFit.a1,0,SiFit.a3,SiFit.a4)).^2;
    Cluster(j).PeakBetasCorrected=(pi/180).*sqrt(Cluster(j).PeakBetasCorrectedSquared); %also convert to radians here
    
    Cluster(j).PeakFWHMsCorrected=(pi/180).*sqrt(Cluster(j).PeakFWHMs.^2-SiFWHM^2);
    Cluster(j).PeakFWHMs2Corrected=(pi/180).*sqrt(Cluster(j).PeakFWHMs2.^2-SiFWHM^2);
    
    %generate Williamson-Hall Plot and fit with line
    LowerBoundWH=[0.0001,lambda/10000];
    UpperBoundWH=[0.01,lambda/50];

    [WHFit,WHGOF]=fit(4*sind(Cluster(j).PeakCenters./2),Cluster(j).PeakBetasCorrected.*cosd(Cluster(j).PeakCenters./2),'poly1','Lower',LowerBoundWH,'Upper',UpperBoundWH);
%     [WHFit2,WHGOF2]=fit(4*sind(Cluster(j).PeakCenters./2),Cluster(j).PeakFWHMsCorrected.*cosd(Cluster(j).PeakCenters./2),'poly1');
%     [WHFit3,WHGOF3]=fit(4*sind(Cluster(j).PeakCenters./2),Cluster(j).PeakFWHMs2Corrected.*cosd(Cluster(j).PeakCenters./2),'poly1');    
    CI=confint(WHFit);
    WHSize = lambda./WHFit.p2; %particle size from WH analysis
    WHStrain = WHFit.p1;
    
%     WHSize2 = lambda./WHFit2.p2;
%     WHStrain2 = WHFit2.p1;
%     
%     WHSize3 = lambda./WHFit3.p2;
%     WHStrain3 = WHFit3.p1;
    
    Cluster(j).WHFit=WHFit;
    Cluster(j).CI=CI;
    Cluster(j).WHGOF=WHGOF;
    Cluster(j).WHSize=WHSize;
    Cluster(j).WHSizeLowCI=lambda./CI(2,2);
    Cluster(j).WHSizeHighCI=lambda./CI(1,2);
    Cluster(j).WHStrain=WHStrain;
    Cluster(j).WHStrainLowCI=CI(1,1);
    Cluster(j).WHStrainHighCI=CI(2,1);
    
%     Cluster(j).WHFit2=WHFit2;
%     Cluster(j).WHGOF2=WHGOF2;
%     Cluster(j).WHSize2=WHSize2;
%     Cluster(j).WHStrain2=WHStrain2;
%     
%     Cluster(j).WHFit3=WHFit3;
%     Cluster(j).WHGOF3=WHGOF3;
%     Cluster(j).WHSize3=WHSize3;
%     Cluster(j).WHStrain3=WHStrain3;
    
    if plot_WH
        figure;
        title(strcat(sampleName,': ', Cluster(j).name,' Williamson-Hall Integral Breadth Plot'))
        hold on
        scatter(4*sind(Cluster(j).PeakCenters./2),Cluster(j).PeakBetasCorrected.*cosd(Cluster(j).PeakCenters./2))
        xlim([0,1])
        ylim([0,0.005])
        text(4*sind(Cluster(j).PeakCenters./2)+0.002,Cluster(j).PeakBetasCorrected.*cosd(Cluster(j).PeakCenters./2)+0.0002,peakLabels2);
        plot(WHFit)
        xlabel('4*Sin(theta)')
        ylabel('Beta*Cos(theta)')
        text(.1,.001,strcat('Particle Size = ',num2str(WHSize),' (Angstroms)','(',num2str(Cluster(j).WHSizeLowCI),',',num2str(Cluster(j).WHSizeHighCI),')'));
        text(.1,.0007,strcat('Particle Strain = ',num2str(WHStrain),'(',num2str(Cluster(j).WHStrainLowCI),',',num2str(Cluster(j).WHStrainHighCI),')'));
        hold off
    
% %         figure;
% %         plot(twoThetaAxis,Cluster(j).DiffProfile)
%         figure;
%         title(strcat(sampleName,' ', Cluster(j).name,' Williamson-Hall FWHM Plot'))
%         hold on
%         scatter(4*sind(Cluster(j).PeakCenters./2),Cluster(j).PeakFWHMsCorrected.*cosd(Cluster(j).PeakCenters./2))
%         xlim([0,1])
%         ylim([0,0.005])
%         text(4*sind(Cluster(j).PeakCenters./2)+0.0005,Cluster(j).PeakFWHMsCorrected.*cosd(Cluster(j).PeakCenters./2)+0.0002,peakLabels);
%         plot(WHFit2)
%         xlabel('4*Sin(theta)')
%         ylabel('Beta*Cos(theta)')
%         text(.6,.001,strcat('Particle Size =',num2str(WHSize2),' ','(Angstroms)'));
%         text(.6,.0007,strcat('Particle Strain =',num2str(WHStrain2)));
%         hold off
%     
% %         figure;
% %         plot(twoThetaAxis,Cluster(j).DiffProfile)
%         figure;
%         title(strcat(sampleName,' ', Cluster(j).name,' Williamson-Hall practical FWHM Plot'))
%         hold on
%         scatter(4*sind(Cluster(j).PeakCenters./2),Cluster(j).PeakFWHMs2Corrected.*cosd(Cluster(j).PeakCenters./2))
%         xlim([0,1])
%         ylim([0,0.005])
%         text(4*sind(Cluster(j).PeakCenters./2)+0.0005,Cluster(j).PeakFWHMs2Corrected.*cosd(Cluster(j).PeakCenters./2)+0.0002,peakLabels);
%         plot(WHFit2)
%         xlabel('4*Sin(theta)')
%         ylabel('Beta*Cos(theta)')
%         text(.6,.001,strcat('Particle Size =',num2str(WHSize3),' ','(Angstroms)'));
%         text(.6,.0007,strcat('Particle Strain =',num2str(WHStrain3)));
%         hold off
    end
end