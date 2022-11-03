% Use this script to process complete diffraction maps of human enamel. The
% script reads in raw 2D diffraction patterns, transforms them into
% cartesian patterns, and analyzes to extract crystallographic parameters
% on a per pattern basis. It outputs .csv files of maps of selected
% parameters for further plotting and analysis with other figure plotting
% scripts

% These analyses support the publication of:
% "Mesoscale Structural Gradients in Human Tooth Enamel", by R. Free, 
% K. DeRocher, V. Cooley, R. Xu, S.R. Stock, and D. Joester.
% This script also includes the units and axes information for each plot.

% Author: Robert Free
% Editor: Derk Joester

%% Clean up
clear all;
close all;
clc
%% Flags
save_maps = true;
save_results = true;
plot_maps = true;

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
c=2.99792*10^18; %angstrom/s
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
sampleName='Enamel D, ';
scanName='50s_pattern';
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

Min2thetaBounds = [14.7639,15.0031,15.2001,15.8234];
Max2thetaBounds = [14.9977,15.1229,15.5055,16.0058];

fit2thetaMinInd=find(twoThetaAxis>14.0,1); %define index for left edge of quadruplet region (used for background subtraction)
fit2thetaMaxInd=find(twoThetaAxis<16.3,1,'last'); %and for right edge
%% Parameter Fitting
%Initialize final plot matrices
AMatrix=zeros(yPts,xPts);
CMatrix=zeros(yPts,xPts);
VolMatrix=zeros(yPts,xPts);
Size121Matrix=zeros(yPts,xPts);
FWHMMatrix=zeros(yPts,xPts);

%Loop over patterns in sample region
for j=pattRange
    %read in background corrected diffractin pattern .tif
    filename=strcat(strcat(scanName,'_',num2str(j,'%04.0f'),'.tif'));
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
    
    %Perform background subtraction around quadruplet
    xDiff=fit2thetaMaxInd-fit2thetaMinInd;
    yDiff=Pattern(j).DiffProfile(fit2thetaMaxInd)-Pattern(j).DiffProfile(fit2thetaMinInd);
    newXaxis=[1:length(twoThetaAxis)]';
    newXaxis=newXaxis-fit2thetaMinInd;
    yAdjustment=(yDiff/xDiff)*newXaxis+Pattern(j).DiffProfile(fit2thetaMinInd);
    AdjustedProfile=Pattern(j).DiffProfile-yAdjustment';
    negativeIndices=find(AdjustedProfile<0); %find indices of all negative values
    AdjustedProfile(negativeIndices)=0; %set negative values to 0
    Pattern(j).DiffProfileQuad=AdjustedProfile; % write to struct
    
    %Perform fit of quadruplet for each center
    
    %initialize storage matrices
    Pattern(j).PeakFitIntensities2=zeros(num2Fit,1);
    Pattern(j).PeakFitCenters2=zeros(num2Fit,1);
    Pattern(j).PeakFitFWHMg2=zeros(num2Fit,1);
    Pattern(j).PeakFitFWHMl2=zeros(num2Fit,1);
    Pattern(j).PeakFitFWHM2=zeros(num2Fit,1);
    Pattern(j).PeakFitEttas=zeros(num2Fit,1);

    %perform fit for pattern in question
    [tempFit,tempGOF]=fit(twoThetaAxis(fit2thetaMinInd:fit2thetaMaxInd)',Pattern(j).DiffProfileQuad(fit2thetaMinInd:fit2thetaMaxInd)',PseudoVoigt4,'StartPoint',StartVals,'Upper',UpperBounds,'Lower',LowerBounds);

    %write fit information to persistent struct
    Pattern(j).fitInfo2=tempFit;
    Pattern(j).GOF2=tempGOF;
    Pattern(j).PeakFitIntensities2=[tempFit.a1,tempFit.a5,tempFit.a9,tempFit.a13];
    Pattern(j).PeakFitCenters2=[tempFit.a2,tempFit.a6,tempFit.a10,tempFit.a14];
    Pattern(j).PeakFitFWHMg2=[tempFit.a3,tempFit.a7,tempFit.a11,tempFit.a15];
    Pattern(j).PeakFitFWHMl2=[tempFit.a4,tempFit.a8,tempFit.a12,tempFit.a16];
    Pattern(j).PeakFitFWHM2=[FWHMcalc(tempFit.a3,tempFit.a4),FWHMcalc(tempFit.a7,tempFit.a8),FWHMcalc(tempFit.a11,tempFit.a12),FWHMcalc(tempFit.a15,tempFit.a16)];
    Pattern(j).PeakFitEttas=[nCalc(tempFit.a3,tempFit.a4),nCalc(tempFit.a7,tempFit.a8),nCalc(tempFit.a11,tempFit.a12),nCalc(tempFit.a15,tempFit.a16)];
    
    %Compute crystallite size from FWHM information
    PeakFitBetas=zeros(num2Fit,1);
    PeakFitBetasCorrected=zeros(num2Fit,1);
    BetaStrain=(180/pi)*4*0.00055*tan((pi/360)*Pattern(j).PeakFitCenters2); %Compute broadening due to average microstrain extracted from W-H of 0.00055. Bstrain = 4*epsilon*tan(theta) 

    %compute integral breadth (area under peak divided by max intensity)
    %for all reflections
    PeakFitBetas=Pattern(j).PeakFitIntensities2./PseudoVoigt(0,Pattern(j).PeakFitIntensities2,0,Pattern(j).PeakFitFWHMg2,Pattern(j).PeakFitFWHMl2);
    smallInd=find(PeakFitBetas<SiFit.a1./PseudoVoigt(0,SiFit.a1,0,SiFit.a3,SiFit.a4)); %determine if Beta is smaller than instrumental value from Si standard
    PeakFitBetasConstrained=PeakFitBetas;
    PeakFitBetasConstrained(smallInd)=(SiFit.a1./PseudoVoigt(0,SiFit.a1,0,SiFit.a3,SiFit.a4))*1.000001; %make sure the value will still be positive after subtraction, but very small.

    %subtract instrumental broadening
    PeakFitBetasCorrectedSquared=PeakFitBetasConstrained.^2-(SiFit.a1./PseudoVoigt(0,SiFit.a1,0,SiFit.a3,SiFit.a4)).^2;
    PeakFitBetasCorrected=sqrt(PeakFitBetasCorrectedSquared);

    %also subtract strain broadening computed through W-H analysis (see
    %below)
    smallInd=find(PeakFitBetasCorrected<BetaStrain);
    PeakFitBetasCorrectedConstrained=PeakFitBetasCorrected;
    PeakFitBetasCorrectedConstrained(smallInd)=1.000001*BetaStrain(smallInd);
    PeakFitBetasCorrectedWH=PeakFitBetasCorrectedConstrained-BetaStrain;

    %Calculate Scherrer using K factors for specific reflections (being
    %careful to convert to radians) and save to persistent struct
    Pattern(j).PeakScherrer=(lambda*Kfactors)./((pi/180)*PeakFitBetasCorrected.*cos((pi/360)*Pattern(j).PeakFitCenters2));
    Pattern(j).PeakScherrerWH=(lambda*Kfactors)./((pi/180)*PeakFitBetasCorrectedWH.*cos((pi/360)*Pattern(j).PeakFitCenters2));

    
    %Compute lattice parameter from mean peak position of 121 and 030
    %reflections

    d003 = lambda/(2*sind(Pattern(j).PeakFitCenters2(3)/2));
    d121 = lambda/(2*sind(Pattern(j).PeakFitCenters2(1)/2));
    Pattern(j).a = 3*d003/cosd(30); %calculate a from 030 reflection
    Pattern(j).c = ((1/d121^2)-((4/3)*7/(Pattern(j).a^2)))^(-1/2); %calculate c from 121 and a parameter
    Pattern(j).Vol=(sqrt(3)/2)*(Pattern(j).a^2)*(Pattern(j).c); %calculate volume for hexagonal crystal system

    
    %Assemble Matrices to visualize across all patterns for each beam
    %center under analysis
    AMatrix(ceil(j/xPts),(1+mod(j-1,xPts)))=Pattern(j).a;
    CMatrix(ceil(j/xPts),(1+mod(j-1,xPts)))=Pattern(j).c;
    VolMatrix(ceil(j/xPts),(1+mod(j-1,xPts)))=Pattern(j).Vol;
    Size121Matrix(ceil(j/xPts),(1+mod(j-1,xPts)))=Pattern(j).PeakScherrerWH(1);
    FWHMMatrix(ceil(j/xPts),(1+mod(j-1,xPts)))=Pattern(j).PeakFitFWHM2(1);
    
    Pattern(j).AzPlots = struct;
    for k=1:num2Fit
        MinBoundInd=find(twoThetaAxis>Min2thetaBounds(k),1);
        MaxBoundInd=find(twoThetaAxis<Max2thetaBounds(k),1,'last');
        Pattern(j).AzPlots(k).data=sum(cart2theta(:,MinBoundInd:MaxBoundInd),2);
        Pattern(j).AzPlots(k).IntIntensity=sum(Pattern(j).AzPlots(k).data);
    end
    
    disp(j)
end
 %% Plot parameter maps

if plot_maps
    fig=figure('Name','a parameter');
    set(fig,'Position',[200,200,500,1000])
    imagesc(AMatrix(:,9:22))
    daspect([yAsp,xAsp,1])
    colormap(gca,parula);
    title('a parameter')
    xlabel(XlabelStr)
    ylabel(YlabelStr)
    cbar = colorbar;
    caxis([9.435,9.448])
    set(gca,'TickDir','out');
    cbar.TickDirection = 'out';

    fig=figure('Name','c parameter');
    set(fig,'Position',[200,200,500,1000])
    imagesc(CMatrix(:,9:22))
    daspect([yAsp,xAsp,1])
    colormap(gca,parula);
    title('c parameter')
    xlabel(XlabelStr)
    ylabel(YlabelStr)
    cbar = colorbar;
    caxis([6.85,6.90])
    set(gca,'TickDir','out');
    cbar.TickDirection = 'out';

    fig=figure('Name','121 Size');
    set(fig,'Position',[200,200,500,1000])
    imagesc(Size121Matrix(:,9:22))
    daspect([yAsp,xAsp,1])
    colormap(gca,parula);
    title('Crystallite Size || 121 (Angstroms)')
    xlabel(XlabelStr)
    ylabel(YlabelStr)
    cbar = colorbar;
    caxis([230,330])
    set(gca,'TickDir','out');
    cbar.TickDirection = 'out';
end
%% Compute CLAA Maps
% Extract Azimuthal profiles from cartesian transforms

% Define 2-theta ranges to bound radial integration for azimuthal plots
% these were defined by reviewing summed diffraction pattern peak edges.
% Could be more systematically extracted from fit.

% Create downsampled patterns for autocorrelation
RedFactor=4; %define factor to downsample at

%rebin the data based on the reduction factor and save in new part of struct
for k=1:num2Fit
    for j=pattRange
        Pattern(j).AzPlots(k).ReducedData=sum(reshape(Pattern(j).AzPlots(k).data,[RedFactor,Nchi/RedFactor]));
    end
end

% Compute autocorrelation
for k=1:num2Fit
    %initialize fields for saving computed data
    AzPlots(k).AutoCorrList=[];
    AzPlots(k).CorrectedAutoCorrList=[];
    AzPlots(k).avgIntList=[];
    AzPlots(k).maxIntList=[];
    
    %compute autocorrelation for each pattern in sample, save important
    %information per plot in AzPlots struct
    for j=pattRange
        [tempAutoCorr,tempLags]=autocorr(Pattern(j).AzPlots(k).ReducedData,'NumLags',((Nchi/RedFactor)-1));
        Pattern(j).AzPlots(k).AutoCorrTotal=tempAutoCorr;
        Pattern(j).AzPlots(k).AutoCorrLags=tempLags;
        Pattern(j).AzPlots(k).AutoCorr=sum(tempAutoCorr(2:9));
        AzPlots(k).AutoCorrList=horzcat(AzPlots(k).AutoCorrList,tempAutoCorr(2));
    end
end

%computed weighted combined azimuthal autocorrelation
for j=pattRange
    
    weightingSum=Pattern(j).AzPlots(1).IntIntensity+Pattern(j).AzPlots(2).IntIntensity+Pattern(j).AzPlots(3).IntIntensity+Pattern(j).AzPlots(4).IntIntensity;
    weights=[Pattern(j).AzPlots(1).IntIntensity/weightingSum,Pattern(j).AzPlots(2).IntIntensity/weightingSum,Pattern(j).AzPlots(3).IntIntensity/weightingSum,Pattern(j).AzPlots(4).IntIntensity/weightingSum];
    AutoCorrs=[Pattern(j).AzPlots(1).AutoCorr,Pattern(j).AzPlots(2).AutoCorr,Pattern(j).AzPlots(3).AutoCorr,Pattern(j).AzPlots(4).AutoCorr];
    Pattern(j).CLAA = dot(weights,AutoCorrs);
    
end

% Arrange autocorrelation values with appropriate x-y relationship
for k=1:num2Fit
    tempCorrMatrix1=zeros(yPts,xPts);

    for j=pattRange
        tempCorrMatrix1(ceil(j/xPts),(1+mod(j-1,xPts)))=Pattern(j).AzPlots(k).AutoCorr;
    end
    AzPlots(k).CorrMatrix=tempCorrMatrix1;
end
CLAAMatrix=zeros(yPts,xPts);
CLAAList=[];
for j=pattRange
    CLAAMatrix(ceil(j/xPts),(1+mod(j-1,xPts)))=Pattern(j).CLAA;
    CLAAList=horzcat(CLAAList,Pattern(j).CLAA);
end

%plot CLAA matrix
if plot_maps
    tempfig=figure;
    imagesc(CLAAMatrix(:,9:22))
    daspect([yAsp,xAsp,1])
    set(tempfig,'Position',[300,400,figWidth,figHeight])
    colormap bone;
    title(strcat(sampleName,'Weighted Combined Autocorrelation, 2D Map'))
    xlabel(XlabelStr)
    ylabel(YlabelStr)
    caxis([0,8])
end
%% Compute average and std deviation of key parameters for ranges of interest

CLAAlist=[];
Size121List=[];
aLatticeList=[];
cLatticeList=[];
rsquaredList=[];
FitParameters=zeros(length(pattRange),17);

for j=pattRange
    CLAAlist=horzcat(CLAAlist,Pattern(j).CLAA);
    Size121List=horzcat(Size121List,Pattern(j).PeakScherrerWH(1));
    aLatticeList=horzcat(aLatticeList,Pattern(j).a);
    cLatticeList=horzcat(cLatticeList,Pattern(j).c);
    rsquaredList=horzcat(rsquaredList,Pattern(j).GOF2.rsquare);
    FitParameters(j,:)=horzcat(j,Pattern(j).fitInfo2.a1,Pattern(j).fitInfo2.a2,Pattern(j).fitInfo2.a3,Pattern(j).fitInfo2.a4,Pattern(j).fitInfo2.a5,Pattern(j).fitInfo2.a6,Pattern(j).fitInfo2.a7,Pattern(j).fitInfo2.a8,Pattern(j).fitInfo2.a9,Pattern(j).fitInfo2.a10,Pattern(j).fitInfo2.a11,Pattern(j).fitInfo2.a12,Pattern(j).fitInfo2.a13,Pattern(j).fitInfo2.a14,Pattern(j).fitInfo2.a15,Pattern(j).fitInfo2.a16);
end

ResultsMatrix = cat(2,pattRange',rsquaredList',CLAAlist',Size121List',aLatticeList',cLatticeList');

CLAAavg=mean(CLAAlist(pattSample));
Size121avg=mean(Size121List(pattSample));
aLatticeavg=mean(aLatticeList(pattSample));
cLatticeavg=mean(cLatticeList(pattSample));


CLAAstd=std(CLAAlist(pattSample));
Size121std=std(Size121List(pattSample));
aLatticestd=std(aLatticeList(pattSample));
cLatticestd=std(cLatticeList(pattSample));
%% Compute the compositional differences in Mg and CO3 corresponding to measured variations in a and c lattice parameters, assuming a simple linear model
ao=9.4554; %from Wilson and colleagues 1999
co=6.8835;

%strain values per substitution based on DeRocher et al 2020
Namg=-0.08;
Ncmg=-0.54;
Naco=-0.22;
Ncco=0.08;

for j=pattSample
    %compute delta a and delta c values relative to sample average
    Pattern(j).DeltaA=Pattern(j).a-aLatticeavg;
    Pattern(j).DeltaC=Pattern(j).c-cLatticeavg;
    Pattern(j).ASoluteExpansion=(Pattern(j).a-aLatticeavg)/aLatticeavg;
    Pattern(j).CSoluteExpansion=(Pattern(j).c-cLatticeavg)/cLatticeavg;    
    
    %compute delta Mg and delta CO3 values relative to sample average
    Pattern(j).DeltaMg=(Pattern(j).DeltaA/(ao*Namg))-(Naco/Namg)*(((Namg*Pattern(j).DeltaC)/(Ncmg*co))-(Pattern(j).DeltaA/ao))/((Namg*Ncco/Ncmg)-Naco);
    Pattern(j).DeltaCarbonate=(((Namg*Pattern(j).DeltaC)/(Ncmg*co))-(Pattern(j).DeltaA/ao))/((Namg*Ncco/Ncmg)-Naco);
   
end

%initialize matrices to hold values
DeltaAMatrix=zeros(yPts,xPts);
DeltaCMatrix=zeros(yPts,xPts);
ASolExpMatrix=zeros(yPts,xPts);
CSolExpMatrix=zeros(yPts,xPts);
DeltaMgMatrix=zeros(yPts,xPts);
DeltaCarbonateMatrix=zeros(yPts,xPts);
%fill matrices for whole sample
for j=pattSample
    DeltaAMatrix(ceil(j/xPts),(1+mod(j-1,xPts)))=Pattern(j).DeltaA;
    DeltaCMatrix(ceil(j/xPts),(1+mod(j-1,xPts)))=Pattern(j).DeltaC;
    ASolExpMatrix(ceil(j/xPts),(1+mod(j-1,xPts)))=Pattern(j).ASoluteExpansion;
    CSolExpMatrix(ceil(j/xPts),(1+mod(j-1,xPts)))=Pattern(j).CSoluteExpansion;
    DeltaMgMatrix(ceil(j/xPts),(1+mod(j-1,xPts)))=Pattern(j).DeltaMg;
    DeltaCarbonateMatrix(ceil(j/xPts),(1+mod(j-1,xPts)))=Pattern(j).DeltaCarbonate;
end
%% plot compositional data
if plot_maps
    tempfig=figure(1);
    imagesc(DeltaAMatrix(:,9:22))
    daspect([yAsp,xAsp,1])
    set(tempfig,'Position',[200,200,figWidth,figHeight])
    title(strcat(sampleName,' Delta A relative to Average'))
    xlabel(XlabelStr)
    ylabel(YlabelStr)
    caxis([-0.005,0.005])

    tempfig=figure(2);
    imagesc(DeltaCMatrix(:,9:22))
    daspect([yAsp,xAsp,1])
    set(tempfig,'Position',[200,200,figWidth,figHeight])
    title(strcat(sampleName,' Delta C relative to Average'))
    xlabel(XlabelStr)
    ylabel(YlabelStr)
    caxis([-0.05,0.05])

    tempfig=figure(3);
    imagesc(DeltaMgMatrix(:,9:22))
    daspect([yAsp,xAsp,1])
    set(tempfig,'Position',[200,200,figWidth,figHeight])
    title(strcat(sampleName,' Delta Mg relative to Average'))
    xlabel(XlabelStr)
    ylabel(YlabelStr)
    caxis([-0.005,0.005])

    tempfig=figure(4);
    imagesc(DeltaCarbonateMatrix(:,9:22))
    daspect([yAsp,xAsp,1])
    set(tempfig,'Position',[200,200,figWidth,figHeight])
    title(strcat(sampleName,' Delta CO3 relative to Average'))
    xlabel(XlabelStr)
    ylabel(YlabelStr)
    caxis([-0.005,0.005])
end
%% Write key results to csv

if save_results
    writematrix(ResultsMatrix,'Fit_Quality_and_Derived_Quantities_3.csv')
    writematrix(FitParameters,'Fit_Parameters_3.csv')
else
    disp('Results Matrix Not Saved!')
end

if save_maps
    %Figure 2
    writematrix(CLAAMatrix,'CLAA_Map.csv')

    %Figure 3
    writematrix(AMatrix,'A_Parameter_Map.csv')
    writematrix(CMatrix,'C_Parameter_Map.csv')
    writematrix(Size121Matrix,'121_Scherrer_Size_Map.csv')

    %Compositional Maps
    writematrix(DeltaAMatrix,'DeltaA_Map.csv')
    writematrix(DeltaCMatrix,'DeltaC_Map.csv')
    writematrix(ASolExpMatrix,'ASoluteExpansion_Map.csv')
    writematrix(CSolExpMatrix,'CSoluteExpansion_Map.csv')
    writematrix(DeltaMgMatrix,'DeltaMg_Map.csv')
    writematrix(DeltaCarbonateMatrix,'DeltaCO3_Map.csv')
else
    disp('Parameter Maps Not Saved!')
end

