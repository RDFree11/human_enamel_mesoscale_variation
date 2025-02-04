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
mfile_name          = 'Fig3_plot.m';

%% Flags
save_figure = true;
save_table = true;

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
c = 2.99792*10^18; %angstrom/s
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
dStep        = (dMax-dMin)/(dRes-1);
twoThetaStep = (twoThetaMax-twoThetaMin)/(twoThetaRes-1);
chiStep      = 360/chiRes;
dAxis        = [dMin:dStep:dMax];
twoThetaAxis = [twoThetaMin:twoThetaStep:twoThetaMax];
chiAxis      = [0:chiStep:360-chiStep];

% Define number of reflections and patterns to analyze
num2Fit       = 4;
numPatt       = 925;
fullPattRange = 1:numPatt;
pattRange     = 1:numPatt;

% Define name of sample and scan (important for reading data)
sampleName = 'Victoria Enamel Sample, ';
scanName   = 'ThinEnamel_fine';

% Define plot options
fig_width     = 9;
fig_asp_ratio = 2;
fig_zoom      = 2;
fig_pos       = [1,1,1+fig_width*fig_zoom,fig_width*fig_zoom/fig_asp_ratio]; %[in]
offset        = 10000;

FontSize_Panel = 10*fig_zoom;
FontSize_Label = 7*fig_zoom;
FontSize_Tick  = 6*fig_zoom;

ang       = char(197); 
cmapstr   = {'bone','winter','winter','winter'};
panl      = {'A','B','C','D','E','F','G','H','I'};
clims     = {[0,8],[27,38],[9.440,9.450],[6.84,6.89]};
dlim_abs  = [[0,8];[27,38];[9.440,9.450];[6.84,6.89]];
cbar_str1 ={{'CLAA','[a.u.]'},{'s_{121}','[nm]'},{'a',['[',ang,']']},{'c',['[',ang,']']}};
cbar_str2 ={{'CLAA*',''},{'s_{121}*','[%]'},{'a*','[%]'},{'c*','[%]'}};
cbar_fmt  = {'','%2.1f','%1.2f','%1.2f'};
cbar_pos  = [0,0.5,0.01,0.3];
rnd       = [1,1,2,2];

%% paths & filenames
pn_csv = './figure source data/';
fn_inp     = {'CLAA_Map',...
              '121_Scherrer_Size_Map',...
              'A_Parameter_Map',...
              'C_Parameter_Map',...
              'Cluster_Map.csv'};

pn_out = './';
fn_out = {'Fig3.pdf','TS2_Crystallographic_Statistics_by_Group.csv'};
          
% Import data
CLAAMatrix=readmatrix('CLAA_Map.csv');
AMatrix=readmatrix('A_Parameter_Map.csv');
CMatrix=readmatrix('C_Parameter_Map.csv');
Size121Matrix=readmatrix('121_Scherrer_Size_Map.csv');

% import cluster map and identify nan patterns
G  = readmatrix(fn_inp{5});
Gp = G(:);
nanPatt = find(isnan(Gp));
GpClean=Gp;
GpClean(nanPatt,:)=[];

% replace Pt patterns with NaN
CLAAMatrix(nanPatt)=NaN;
AMatrix(nanPatt)=NaN;
CMatrix(nanPatt)=NaN;
Size121Matrix(nanPatt)=NaN;

%assemble into single matrix
M(:,:,1) = CLAAMatrix;
M(:,:,2) = Size121Matrix;
M(:,:,3) = AMatrix;
M(:,:,4) = CMatrix;

% select data and convert units
M = M(:,x_range,:);
M(:,:,2) = M(:,:,2)/10; % convert from Angstrom to nm

% convert to 2D array
Mp = squeeze(reshape(M,[],1,length(fn_inp)-1)); % reshape to vector

% remove NaN values to compute sample average
MpClean = Mp;
MpClean(nanPatt,:)=[];

% calculate overall mean and standard deviation
mean_val = mean(MpClean);
std_val  = std(MpClean);

Mp(nanPatt,:)=repmat([NaN,NaN,NaN,NaN],length(nanPatt),1);

dlim_rel  = 100*(dlim_abs./mean_val'-1);

%calculate relative deviation from mean value
for ii = 1:4
    C(:,:,ii) = M(:,:,ii)/mean_val(ii)-1;
end

%% compute statistics by group
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

%% compute statistics for whole sample
sample_mean = mean(MpClean,1);
sample_std  = std(MpClean,1);
sample_median = median(MpClean,1);
sample_quantiles = quantile(MpClean,[.25,.75],1);
sample_iqr = sample_quantiles(2,:)-sample_quantiles(1,:);
sample_lower_median = sample_median - 1.57.*(sample_iqr)./sqrt(numel(GpClean));
sample_upper_median = sample_median + 1.57.*(sample_iqr)./sqrt(numel(GpClean));

outlierIndices=zeros(numel(GpClean),4); %initialize binary array to track outlier indices
for ii=1:4
    for jj=1:numel(GpClean)
        if (MpClean(jj,ii)>sample_quantiles(2,ii)+1.5*sample_iqr(ii))
            outlierIndices(jj,ii)=1;
        end
        if (MpClean(jj,ii)<sample_quantiles(1,ii)-1.5*sample_iqr(ii))
            outlierIndices(jj,ii)=1;
        end
    end
end

% replace outlier values with NaN
Mp_noOutliers_sample=MpClean;
Mp_noOutliers_sample(find(outlierIndices==1))=NaN;

%calculate sample max and min without outlier values removed
sample_max=max(Mp_noOutliers_sample,[],1);
sample_min=min(Mp_noOutliers_sample,[],1);

%% Plot
close all
xv = [0,diff(x_range([1,end]))*dx];
yv = [0,(size(M,1)-1)*dy];

fig=figure;
set(fig,'Units','inches','Position',fig_pos)


cbarxoffset = 0.025;
cbaryoffset = 0.68;

for ii=1:4
    ax(ii) = subplot(2,30,[(ii-1)*6+[1:4]]);
    imagesc(xv,yv,M(:,:,ii));
    text(0,1.15,panl{ii},'Units','normalized','FontSize',FontSize_Panel,'Color','k');
    ax(ii).TickDir = 'out';
    colormap(ax(ii),cmapstr{ii});
    cbar(ii) = colorbar;
    cbar(ii).FontSize = FontSize_Tick;
    cbar(ii).TickDirection = 'out';
    ap = ax(ii).Position;
    if ii>1    
        cbar(ii).Position = [ap(1)+ap(3)+cbarxoffset,cbaryoffset,0.005,0.15];
        text(ax(ii),1.4,1.25,cbar_str1{ii},...
            'Units','normalized',...
            'FontSize',FontSize_Tick,...
            'Color','k',...
            'HorizontalAlignment','center');
        text(ax(ii),1.15,1.25,cbar_str2{ii},...
            'Units','normalized',...
            'FontSize',FontSize_Tick,...
            'Color','k',...
            'HorizontalAlignment','center');
    else
        cbar(ii).Position = [ap(1)+ap(3)+cbarxoffset-0.015,cbaryoffset,0.005,0.15];
        text(ax(ii),1.2,1.25,cbar_str1{ii},...
            'Units','normalized',...
            'FontSize',FontSize_Tick,...
            'Color','k',...
            'HorizontalAlignment','center');
    end
       
    daspect([yAsp,xAsp,1])
    ax(ii).FontSize = FontSize_Tick;
    xlabel(XlabelStr,'FontSize',FontSize_Label);
    ax(ii).TickDir = 'out';
    ax(ii).CLim = clims{ii};
end

% manually label clusters. Note: this may change as k-means 
% clustering assigns cluster number 
% cats = categorical(["RH","RT","IR"]);
% cats = reordercats(cats,["RH","RT","IR"]);

cats = categorical(["RH","IR"]);
cats = reordercats(cats,["RH","IR"]);

% plot CLAA map
ax(5) = subplot(2,30,[25:28]);
imagesc(xv,yv,M(:,:,1));
ax(5).TickDir  = 'out';
ax(5).YTick    = [];
ax(5).FontSize = FontSize_Tick;
daspect([yAsp,xAsp,1])
colormap(ax(5),'bone');
xlabel(XlabelStr,'FontSize',FontSize_Label);


ax(6)=axes('Parent',fig,...
           'Units',ax(5).Units,'Position',ax(5).Position,...
           'Color','none');

% determine transparency map 
tmap = ones(size(G));
% set colormap
% cmap = jet(3);
cmap = [0 0 1;0 .75 0];

% superpose transparent cluster map
alpha_val = 0.25;
imagesc(ax(6),xv,yv,G,'alphadata',alpha_val*tmap);
daspect([yAsp,xAsp,1])
colormap(ax(6),cmap);
caxis(ax(6),[0.5,3.5]);
ax(6).Visible = 'off';
linkprop([ax(5),ax(6)],'Position');

% draw and format colorbar
cbar(1) = colorbar;
% cbar(1).Ticks  = 1:3;
% cbar(1).TickLabels = {'RH','RT','IR'};
cbar(1).Ticks  = 1:2:3;
cbar(1).TickLabels = {'RH','IR'};
cbar(1).Direction = 'reverse';
cbar(1).TickDirection = 'out';
cbar(1).FontSize = FontSize_Tick;
ap = ax(6).Position;
cbar(1).Position = [ap(1)+ap(3)+cbarxoffset-0.015,cbaryoffset,0.005,0.15];
text(0,1.15,panl{5},'Units','normalized','FontSize',FontSize_Panel,'Color','k');

% add second axis to colorbars for B-D
for ii=2:4
    ax(ii).YTick = [];
    yt = round(100*diff(clims{ii})/mean_val(ii)/5,rnd(ii))*[-2:1:2];
    yl = 100*(clims{ii}/mean_val(ii)-1); % [%]
    ytlabel = num2str(yt',cbar_fmt{ii});
    new_ax(ii)=axes(...
        'Parent',fig,...
        'Units',cbar(ii).Units,'Position',cbar(ii).Position,...
        'Color','none',...
        'YTick',yt,...
        'YLim', yl,...
        'YTickLabels',ytlabel,...
        'XTick',[],...
        'TickDir','out',...
        'FontSize',FontSize_Tick);
    end

%% boxplots for panels F-I
for ii=1:4
    if ii==1
        ax(ii+5) = subplot(2,30,8*(ii-1)+30+[1:6]);
    end
    if ii==2
        ax(ii+5) = subplot(2,30,8*(ii-1)+29+[1:7]);        
    end
    if ii==3
        ax(ii+5) = subplot(2,30,8*(ii-1)+30+[1:7]);        
    end
    if ii==4
        ax(ii+5) = subplot(2,30,8*(ii-1)+30+[1:6]);        
    end
    boxplot(MpClean(:,ii),cats(GpClean),...
            'Notch','on',...
            'DataLim',dlim_abs(ii,:),...
            'ExtremeMode','clip');hold on; %

    % fill boxes using colors from group map
    % note that the order in h is the reverse of the axis for some reason
    h = findobj(ax(ii+5),'Tag','Box');
    for kk=1:length(h)
         patch(get(h(kk),'XData'),get(h(kk),'YData'),cmap(end+1-kk,:));
    end

    % annotate with group mean
    plot(group_mean(:,ii),'r_')
    % format plot
    ax(ii+5).YMinorTick = 'on';
    ax(ii+5).TickDir = 'out';
    ax(ii+5).FontSize = FontSize_Tick;

    % indicate mean and standard deviation across all clusters
    yline(mean_val(ii),'k-','LineWidth',1);
    yline(mean_val(ii)+std_val(ii),'k:','LineWidth',1);
    yline(mean_val(ii)-std_val(ii),'k:','LineWidth',1);
    xlabel('cluster ID','FontSize',FontSize_Label)
    ylabel(join(cbar_str1{ii}),'FontSize',FontSize_Label);
    text(ax(ii+5),0,1.1,panl{ii+5},'FontSize',FontSize_Panel,'Color','k','Units','normalized')
    yyaxis right
    ax(ii+5).YAxis(2).Limits=100*(ax(ii+5).YAxis(1).Limits/mean_val(ii)-1);
    ax(ii+5).YAxis(2).MinorTick = 'on';
    ylabel(join(cbar_str2{ii}),'FontSize',FontSize_Label)
end

%% Tabulate Results By Cluster
clc
TS2=table('Size',[0,11],...
          'VariableTypes',{'categorical','categorical','double','double','double','double','double','double','double','double','double'},...
          'VariableNames',{'Varname','Cluster Group','Mean','StdDev','Non-Outlier Max','Upper Quartile',...
          '95% CI of Medial Upper Bound','Median','95% CI of Medial Lower Bound','Lower Quartile','Non-Outlier Minimum'});

varNames = {'CLAA','s121','a Parameter','c Parameter'};
% clusterName = {['RH n=',num2str(group_nnz(1,1))],['RT n=',num2str(group_nnz(2,1))],['IR n=',num2str(group_nnz(3,1))],'Sample'};

%build table from group statistics for 3 clusters
% for ii=1:4 %loop over variables
%     for jj=1:4 %loop over groups and sample
%         if jj<4
%             TS2(4*(ii-1)+jj,:)={varNames{ii},clusterName{jj},group_mean(jj,ii),group_std(jj,ii),...
%                 group_max(jj,ii),group_quantiles(jj,ii,2),group_upper_median(jj,ii),group_median(jj,ii),...
%                 group_lower_median(jj,ii),group_quantiles(jj,ii,1),group_min(jj,ii)};
%         else
%             TS2(4*(ii-1)+jj,:)={varNames{ii},clusterName{jj},sample_mean(ii),sample_std(ii),...
%                 sample_max(ii),sample_quantiles(2,ii),sample_upper_median(ii),sample_median(ii),...
%                 sample_lower_median(ii),sample_quantiles(1,ii),sample_min(ii)};
%         end        
%     end    
% end


clusterName = {['RH n=',num2str(group_nnz(1,1))],['IR n=',num2str(group_nnz(2,1))],'Sample'};
%build table from group statistics for 2 clusters
for ii=1:4 %loop over variables
    for jj=1:3 %loop over groups and sample
        if jj<3
            TS2(3*(ii-1)+jj,:)={varNames{ii},clusterName{jj},group_mean(jj,ii),group_std(jj,ii),...
                group_max(jj,ii),group_quantiles(jj,ii,2),group_upper_median(jj,ii),group_median(jj,ii),...
                group_lower_median(jj,ii),group_quantiles(jj,ii,1),group_min(jj,ii)};
        else
            TS2(3*(ii-1)+jj,:)={varNames{ii},clusterName{jj},sample_mean(ii),sample_std(ii),...
                sample_max(ii),sample_quantiles(2,ii),sample_upper_median(ii),sample_median(ii),...
                sample_lower_median(ii),sample_quantiles(1,ii),sample_min(ii)};
        end        
    end    
end

%% Save Results

if save_figure
    saveas(gcf,[pn_out,fn_out{1}],'pdf');
    disp(['Saved as: ',pn_out,fn_out{1}])
else
    disp('Figure not saved!');
end

if save_table
    writetable(TS2,[pn_out,fn_out{2}])
    disp(['Saved as: ',pn_out,fn_out{2}])
else
    disp('Table not saved!')
end
    