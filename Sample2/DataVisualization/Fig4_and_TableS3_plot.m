% Use this script in conjunction with the .csv data files to generate rough
% versions of the figures and tables for 
% "Mesoscale Structure and Composition Varies Systematically in Human Tooth
% Enamel", by R. Free, K. DeRocher, V. Cooley, R. Xu, S.R. Stock, and D. Joester.
% This script also includes the units and axes information for each plot.

% Author: Robert Free
% Editor: Derk Joester

%% Clean Up
clear variables 
close all
clc

%% Set path (run script rather than section) 
% set path
mfile_name          = 'Fig4_plot.m';

%% Flags
save_figure = true;
save_tables = true;

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
x_range = 1:32;

% labels for CLAA map
XlabelStr = 'x [µm]';
YlabelStr = 'y [µm]';

% Define dimensions of map in pixels (different for each scan)
xPts = 37;
yPts = 25;

energyKeV = 17.00000; %energy for scan

%physical constants
h = 4.135667*10^-18; % keV*s
c = 2.998792*10^18;  % angstrom/s
lambda = h*c/energyKeV; % angstrom

samp2det = 155.6166; %sample to detector distance in mm (may change per scan)
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
dStep=(dMax-dMin)/(dRes-1);
twoThetaStep=(twoThetaMax-twoThetaMin)/(twoThetaRes-1);
chiStep=360/chiRes;
dAxis=dMin:dStep:dMax;
twoThetaAxis=twoThetaMin:twoThetaStep:twoThetaMax;
chiAxis=0:chiStep:360-chiStep;

% Define number of reflections and patterns to analyze
num2Fit=4;
numPatt=925;
fullPattRange=1:numPatt;
pattRange=1:numPatt;

% Define name of sample and scan (important for reading data)
sampleName = 'Victoria Enamel Sample, ';
scanName   = 'ThinEnamel_fine';

% Define plot options
fig_width     = 10;
fig_asp_ratio = 2;
fig_zoom      = 2;
fig_pos = [1,1,1+fig_width*fig_zoom,fig_width*fig_zoom/fig_asp_ratio]; %[in]
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
pn_csv = './';
fn_inp     = {'CLAA_Map',...
              'DeltaC_Map',...
              'DeltaA_Map',...
              'DeltaMg_Map'...
              'DeltaCO3_Map'};

pn_out = './';
fn_out = {'Fig4.eps','Fig4_Table_Means.csv','Fig4_Table_MultComparison.csv'};
          
% Import data
CLAAMatrix=readmatrix('CLAA_Map.csv');
DeltaCMatrix=readmatrix('DeltaC_Map.csv');
DeltaAMatrix=readmatrix('DeltaA_Map.csv');
DeltaMgMatrix=readmatrix('DeltaMg_Map.csv');
DeltaCO3Matrix=readmatrix('DeltaCO3_Map.csv');

% import cluster map and identify nan patterns
G  = readmatrix('Cluster_Map.csv');
Gp = G(:);
nanPatt = find(isnan(Gp));
GpClean=Gp;
GpClean(nanPatt,:)=[];

% replace Pt patterns with NaN
CLAAMatrix(nanPatt)=NaN;
DeltaCMatrix(nanPatt)=NaN;
DeltaAMatrix(nanPatt)=NaN;
DeltaMgMatrix(nanPatt)=NaN;
DeltaCO3Matrix(nanPatt)=NaN;


M(:,:,1) = CLAAMatrix;
M(:,:,2) = DeltaCMatrix;
M(:,:,3) = DeltaAMatrix;
M(:,:,4) = DeltaMgMatrix;
M(:,:,5) = DeltaCO3Matrix;

% convert from fraction to %
M(:,:,2:5) = 100*M(:,:,2:5);
% crop to ROI
M = M(:,x_range,:);

%convert to 2D array
Mp = squeeze(reshape(M,[],1,length(fn_inp))); % reshape to vector

% remove NaN values to compute sample average
nanPatt = find(isnan(Mp(:,1)));
MpClean = Mp;
MpClean(nanPatt,:)=[];

% calculate overall mean and standard deviation
mean_val = mean(MpClean);
std_val  = std(MpClean);

Mp(nanPatt,:)=repmat([NaN,NaN,NaN,NaN,NaN],length(nanPatt),1);

% axis vectors
xv = [0,diff(x_range([1,end]))*dx];
yv = [0,(size(M(:,:,1),1)-1)*dy];

%% Compute Statistics by group
[group_nnz,group_mean,group_std,group_median,group_quantiles]=grpstats(MpClean,GpClean,{'nnz','mean','std','median',@(x)quantile(x,[.25,.75])});

% compute 95% confidence interval from median as boxplot function does
group_iqr = group_quantiles(:,:,2) - group_quantiles(:,:,1);
group_lower_median = group_median - 1.57.*(group_iqr)./sqrt(group_nnz);
group_upper_median = group_median + 1.57.*(group_iqr)./sqrt(group_nnz);

% deterimine outliers as points more than 1.5x the interquartile range
% below or above the 25th or 75th percentiles, respectively.
outlierIndices=zeros(numel(GpClean),4); %initialize binary array to track outlier indices
for ii=1:4
    for jj=1:numel(GpClean)
        if (MpClean(jj,ii)>group_quantiles(GpClean(jj),ii,2)+1.5*group_iqr(GpClean(jj),ii))
            outlierIndices(jj,ii)=1;
        end
        if (MpClean(jj,ii)<group_quantiles(GpClean(jj),ii,1)-1.5*group_iqr(GpClean(jj),ii))
            outlierIndices(jj,ii)=1;
        end
    end
end

% replace outlier values with NaN
Mp_noOutliers_group=MpClean;
Mp_noOutliers_group(find(outlierIndices==1))=NaN;

%calculate max and min without outlier values removed
[group_max,group_min]=grpstats(Mp_noOutliers_group,GpClean,{'max','min'});

group_sem   = group_std./group_nnz;
group_mean_rel = group_mean./mean_val-1;
group_std_rel  = group_std./mean_val;
group_sem_rel  = group_sem./mean_val;
%% Plot
close all
% manually label clusters. Note: this may change as k-means 
% clustering assigns cluster number 
% cats = categorical(["RH","RT","IR"]);
% cats = reordercats(cats,["RH","RT","IR"]);

cats = categorical(["RH","IR"]);
cats = reordercats(cats,["RH","IR"]);

cbarxoffset = 0.005;
cbaryoffset = 0.355;

fig = figure('Unit','inches','Position',fig_pos);
% plot CLAA map
ax(1) = subplot(2,16,[1:3,17:19]);
imagesc(xv,yv,M(:,:,1));
ax(1).TickDir = 'out';
daspect([yAsp,xAsp,1])
colormap(ax(1),'bone');
xlabel(XlabelStr,'FontSize',FontSize_Label);
ylabel(YlabelStr,'FontSize',FontSize_Label);

ax(2)=axes('Parent',fig,...
           'Units',ax(1).Units,'Position',ax(1).Position,...
           'Color','none');

% determine transparency map 
tmap = ones(size(G));
% set colormap
% cmap = jet(3);
cmap = [0 0 1;1 1 0];

% superpose transparent cluster map
imagesc(ax(2),xv,yv,G,'alphadata',alpha_val*tmap);
daspect([yAsp,xAsp,1])
colormap(ax(2),cmap);
caxis(ax(2),[0.5,3.5]);
ax(2).Visible = 'off';
linkprop([ax(1),ax(2)],'Position');

% draw and format colorbar
cbar(1) = colorbar;
% cbar(1).Ticks  = 1:3;
% cbar(1).TickLabels = {'RH','RT','IR'};
cbar(1).Ticks  = 1:2:3;
cbar(1).TickLabels = {'RH','IR'};
cbar(1).Direction = 'reverse';
cbar(1).TickDirection = 'out';
ap = ax(2).Position;
cbar(1).Position = [ap(1)+ap(3)+cbarxoffset,.455,0.01,0.1];
text(ax(2),0,1.3,panl{1},'FontSize',FontSize_Panel,'Color','k','Units','normalized')

% plot Mg map and Carbonate Map only
for ii = 2:3
    ax(ii) = subplot(2,16,4*(ii-1)+[1:3,17:19]);
    imagesc(xv,yv,M(:,:,ii+2)); 
    ax(ii).TickDir = 'out';
    daspect([yAsp,xAsp,1])
    colormap(ax(ii),'parula');
    xlabel(XlabelStr,'FontSize',FontSize_Label);
    ylabel(YlabelStr,'FontSize',FontSize_Label);
    caxis(ax(ii),[-0.6,0.6]);
    cbar(ii) = colorbar;
    cbar(ii).TickDirection = 'out';
    ap = ax(ii).Position;
    cbar(ii).Position = [ap(1)+ap(3)+cbarxoffset,cbaryoffset,0.01,0.3];
    text(ax(ii),.9,1.35,cbar_str{ii},...
    'Units','normalized',...
    'FontSize',FontSize_Tick,...
    'Color','k',...
    'HorizontalAlignment','center');

    text(ax(ii),0,1.3,panl{ii},'FontSize',FontSize_Panel,'Color','k','Units','normalized')
end

%% plot group means for predicted Mg map as a bar plot
ii = 4;
ax(ii) = subplot(2,16,[13:16]);
boxplot(MpClean(:,ii),cats(GpClean),...
        'Notch','on',...
        'DataLim',[-0.6,0.6],...
        'ExtremeMode','clip');hold on; %

% fill boxes using colors from group map
% note that the order in h is the reverse of the axis for some reason
h = findobj(ax(ii),'Tag','Box');
for kk=1:length(h)
     patch(get(h(kk),'XData'),get(h(kk),'YData'),cmap(end+1-kk,:));
end

% annotate with group mean
plot(group_mean(:,ii),'r_')
% format plot
ax(ii).YMinorTick = 'on';
ax(ii).TickDir = 'out';
ax(ii).FontSize = FontSize_Tick;

% indicate mean and standard deviation across all clusters
yline(mean_val(ii),'k-','LineWidth',1);
yline(mean_val(ii)+std_val(ii),'k:','LineWidth',1);
yline(mean_val(ii)-std_val(ii),'k:','LineWidth',1);
xlabel('cluster ID','FontSize',FontSize_Label)
text(ax(ii),0,1.1,panl{ii},'FontSize',FontSize_Panel,'Color','k','Units','normalized')

% plot group means for predicted carbonate map as a bar plot
ax(ii+1) = subplot(2,16,[29:32]);

boxplot(MpClean(:,ii+1),cats(GpClean),...
        'Notch','on',...
        'DataLim',[-0.6,0.6],...
        'ExtremeMode','clip');hold on; %
    
% fill boxes using colors from group map
% note that the order in h is the reverse of the axis for some reason
h = findobj(ax(ii+1),'Tag','Box');
for kk=1:length(h)
     patch(get(h(kk),'XData'),get(h(kk),'YData'),cmap(end+1-kk,:));
end
% annotate with mean over all clusters
plot(group_mean(:,ii+1),'r_')
ax(ii+1).TickDir = 'out';
ax(ii+1).YMinorTick = 'on';
ax(ii+1).FontSize = FontSize_Tick;

% indicate standard deviation across all clusters
yline(0,'k-','LineWidth',1);
yline(std_val(ii+1),'k:','LineWidth',1);
yline(-std_val(ii+1),'k:','LineWidth',1); 
xlabel('cluster ID','FontSize',FontSize_Label)
text(ax(ii+1),0,1.1,panl{ii+1},'FontSize',FontSize_Panel,'Color','k','Units','normalized')

if save_figure
    saveas(gcf,[pn_out,fn_out{1}],'pdf');
    disp(['Saved as: ',pn_out,fn_out{1}])
else
    disp('Figure not saved!');
end
    
%% 1-way ANOVA & Multiple Comparison
clc
T1 = table('Size',[0,7],...
          'VariableTypes',{'categorical','categorical','categorical','double','double','double','double'},...
          'VariableNames',{'Varname','Group1','Group2','LowerCI','Difference','UpperCI','p'});
      
for ii = 1:size(Mp,2)
    [~,~,stats] = anova1(MpClean(:,ii),cats(GpClean),'off');
    [c,m,~,gnames] = multcompare(stats,'Display','off');
    Temp = array2table(c,'VariableNames',{'Group1','Group2','LowerCI','Difference','UpperCI','p'});
    Temp.Varname = categorical(cellstr(repmat(varname{ii},height(Temp),1)));
    Temp.Group1 = categorical(gnames(Temp.Group1));
    Temp.Group2 = categorical(gnames(Temp.Group2));
    T1 = [T1;Temp];
end
disp(T1)

T2=table('Size',[0,11],...
          'VariableTypes',{'categorical','categorical','double','double','double','double','double','double','double','double','double'},...
          'VariableNames',{'Varname','Cluster Group','Mean','StdDev','Non-Outlier Max','Upper Quartile',...
          '95% CI of Medial Upper Bound','Median','95% CI of Medial Lower Bound','Lower Quartile','Non-Outlier Minimum'});

varNames = {'DeltaC/C','DeltaA/A','DeltaMg','DeltaCO3'};

clusterName = {['RH n=',num2str(group_nnz(1,1))],['IR n=',num2str(group_nnz(2,1))],'Sample'};

%build table from group statistics
for ii=1:4 %loop over variables
    for jj=1:2 %loop over groups and sample
        T2(2*(ii-1)+jj,:)={varNames{ii},clusterName{jj},group_mean(jj,ii+1),group_std(jj,ii+1),...
            group_max(jj,ii+1),group_quantiles(jj,ii+1,2),group_upper_median(jj,ii+1),group_median(jj,ii+1),...
            group_lower_median(jj,ii+1),group_quantiles(jj,ii+1,1),group_min(jj,ii+1)};  
    end    
end
disp(T2)

if save_tables
    writetable(T1,[pn_out,fn_out{3}]);
    writetable(T2,[pn_out,fn_out{2}]);
    disp(['Saved as: ',pn_out,fn_out{2}])
else
    disp('Figure not saved!');
end