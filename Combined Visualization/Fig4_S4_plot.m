% Use this script in conjunction with the .csv data files to generate rough
% versions of the figures for 
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

% define pixel pitch for CLAA maps
dx = 0.5; % [µm]
dy = 0.5; % [µm]

% crop CLAA map to indices for each sample
Sample(1).x_range = 4:52;
Sample(2).x_range = 1:32;
Sample(3).x_range = 9:22;

% labels for CLAA map
XlabelStr = 'x [µm]';
YlabelStr = 'y [µm]';

% Define dimensions of map in pixels for each scan
Sample(1).xPts = 52;
Sample(1).yPts = 15;
Sample(2).xPts = 37;
Sample(2).yPts = 25;
Sample(3).xPts = 23;
Sample(3).yPts = 41;

energyKeV = 17.00000; %energy for scans

%physical constants
h = 4.135667*10^-18; % keV*s
c = 2.998792*10^18;  % angstrom/s
lambda = h*c/energyKeV; % angstrom

Sample(1).samp2det  = 155.3818; % sample to detector distance in mm (may change per scan)
Sample(2).samp2det = 155.6166;
Sample(3).samp2det = 109.2134;

pixelSize = 0.079278; % pixel size in mm (usually the same for MAR detector from 34ide)

% Define name of sample and scan (important for reading data)
Sample(1).name="Sample1";
Sample(2).name="Sample2";
Sample(3).name="Sample3";

% Define plot options
fig_width     = 12;
fig_asp_ratio = 1;
fig_zoom      = 1;
fig_pos       = [1,1,1+fig_width*fig_zoom,fig_width*fig_zoom/fig_asp_ratio]; %[in]
offset        = 10000;

FontSize_Panel = 10*fig_zoom;
FontSize_Label = 7*fig_zoom;
FontSize_Tick  = 6*fig_zoom;

ang  = char(197); 
panl = {'A','B','C','D','E'};
varname = {'\Delta_a/a','\Delta_c/c','\Deltax_Mg','\Deltax_CO3'};

yl ={4+[-2,2],29+[-2.5,2.5],9.447+[-0.0025,0.0025],6.880+[-0.015,0.015]};


cbar_str = {'\Deltaa / a_{avg} [%]',...
              {'\Deltac / c_{avg} [%]'},...
              '\Deltax_{Mg} [at%]',...
              '\Deltax_{CO_3} [at%]'};
alpha_val = 0.25;

ang       = char(197); 
cmapstr   = {'autumn','winter','summer';'autumn','winter','summer';'autumn','winter','summer';'autumn','winter','summer'};
% cmapstr   = {'bone','bone','bone';'autumn','autumn','autumn';'autumn','autumn','autumn';'autumn','autumn','autumn'};
panl      = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S'};
climsGlobal = {[-0.060,0.060],[-.40,.40],[-0.8,0.8],[-0.5,0.5],[0,8]};
climsSample = {[-0.060,0.060],[-.40,.40],[-0.8,0.8],[-0.5,0.5],[0,8];[-0.055,0.055],[-0.40,0.40],[-0.5,0.5],[-0.5,0.5],[0,8];[-0.055,0.055],[-0.35,0.35],[-0.5,0.5],[-0.5,0.5],[0,8]};
dlim_abs  = [[-0.060,0.060];[-.40,.40];[-0.8,0.8];[-0.5,0.5];[0,8]];
cbar_fmt  = {'','%2.1f','%1.2f','%1.2f'};
cbar_pos  = [0,0.5,0.01,0.3];
rnd       = [1,1,2,2];

%% Read in data

% paths & filenames
fn_inp     = {'ASoluteExpansion_Map',...
              'CSoluteExpansion_Map',...
              'DeltaMg_Map'...
              'DeltaCO3_Map'...
              'CLAA_Map'};

pn_out = './';
fn_out = {'Fig4.pdf','Fig4_Table_Means.csv','Fig4_Table_MultComparison.csv'};
          
% Import data
for kk=1:3
    
    DeltaAMatrix=zeros(Sample(kk).yPts,Sample(kk).xPts);
    DeltaCMatrix=zeros(Sample(kk).yPts,Sample(kk).xPts);
    DeltaMgMatrix=zeros(Sample(kk).yPts,Sample(kk).xPts);
    DeltaCO3Matrix=zeros(Sample(kk).yPts,Sample(kk).xPts);
    CLAAMatrix=zeros(Sample(kk).yPts,Sample(kk).xPts);
    
    folderString = Sample(kk).name;
    DeltaAMatrix=readmatrix(strcat(folderString,'/','ASoluteExpansion_Map.csv'));
    DeltaCMatrix=readmatrix(strcat(folderString,'/','CSoluteExpansion_Map.csv'));
    DeltaMgMatrix=readmatrix(strcat(folderString,'/','DeltaMg_Map.csv'));
    DeltaCO3Matrix=readmatrix(strcat(folderString,'/','DeltaCO3_Map.csv'));
    CLAAMatrix=readmatrix(strcat(folderString,'/','CLAA_Map'));
    
    % import cluster map and identify NaN patterns
    G=readmatrix(strcat(folderString,'/','Cluster_Map.csv'));
    Gp=G(:);
    nanPatt = find(isnan(Gp));
    GpClean=Gp;
    GpClean(nanPatt,:)=[];
    
    % replace Pt patterns with NaN
    
    DeltaAMatrix(nanPatt)=NaN;
    DeltaCMatrix(nanPatt)=NaN;
    DeltaMgMatrix(nanPatt)=NaN;
    DeltaCO3Matrix(nanPatt)=NaN;
    CLAAMatrix(nanPatt)=NaN;
    
    %assemble into single Matrix
    clear M
    M(:,:,1) = DeltaAMatrix;
    M(:,:,2) = DeltaCMatrix;
    M(:,:,3) = DeltaMgMatrix;
    M(:,:,4) = DeltaCO3Matrix;
    M(:,:,5) = CLAAMatrix;
    
    % select data and convert from fraction to %
    M = M(:,Sample(kk).x_range,:);
    M(:,:,1:4) = 100*M(:,:,1:4);
    
    %convert to 2D array
    Mp = squeeze(reshape(M,[],1,length(fn_inp))); % reshape to vector
    
    % remove NaN values to compute sample average
    MpClean = Mp;
    MpClean(nanPatt,:)=[];
    Mp(nanPatt,:)=repmat([NaN,NaN,NaN,NaN,NaN],length(nanPatt),1);

    % save key variables to sample struct   
    Sample(kk).M=M;
    Sample(kk).Mp=Mp;
    Sample(kk).MpClean=MpClean;
    Sample(kk).G=G;
    Sample(kk).GpClean=GpClean;

    % calculate sample mean and standard deviation
    Sample(kk).mean_val=mean(MpClean);
    Sample(kk).std_val=std(MpClean);
    
end
%% compute statistics by group for each sample

for kk=1:3

    [group_nnz,group_mean,group_std,group_median,group_quantiles]=grpstats(Sample(kk).MpClean,Sample(kk).GpClean,{'nnz','mean','std','median',@(x)quantile(x,[.25,.75])});
    
    % compute 95% confidence interval from median as boxplot function does
    group_iqr = group_quantiles(:,:,2) - group_quantiles(:,:,1);
    group_lower_median = group_median - 1.57.*(group_iqr)./sqrt(group_nnz);
    group_upper_median = group_median + 1.57.*(group_iqr)./sqrt(group_nnz);
    
    % deterimine outliers as points more than 1.5x the interquartile range
    % below or above the 25th or 75th percentiles, respectively.
    outlierIndices=zeros(numel(Sample(kk).GpClean),5); %initialize binary array to track outlier indices
    for ii=1:5
        for jj=1:numel(Sample(kk).GpClean)
            if (Sample(kk).MpClean(jj,ii)>group_quantiles(Sample(kk).GpClean(jj),ii,2)+1.5*group_iqr(Sample(kk).GpClean(jj),ii))
                outlierIndices(jj,ii)=1;
            end
            if (Sample(kk).MpClean(jj,ii)<group_quantiles(Sample(kk).GpClean(jj),ii,1)-1.5*group_iqr(Sample(kk).GpClean(jj),ii))
                outlierIndices(jj,ii)=1;
            end
        end
    end
    
    % replace outlier values with NaN
    Mp_noOutliers_group=Sample(kk).MpClean;
    Mp_noOutliers_group(find(outlierIndices==1))=NaN;
    
    %calculate max and min without outlier values removed
    [group_max,group_min]=grpstats(Mp_noOutliers_group,Sample(kk).GpClean,{'max','min'});
    
    group_sem   = group_std./group_nnz;
    group_mean_rel = group_mean./Sample(kk).mean_val-1;
    group_std_rel  = group_std./Sample(kk).mean_val;
    group_sem_rel  = group_sem./Sample(kk).mean_val;

    % Write key variables to sample struct
    Sample(kk).group_nnz=group_nnz;
    Sample(kk).group_mean=group_mean;
    Sample(kk).group_std=group_std;
    Sample(kk).group_median=group_median;
    Sample(kk).group_quantiles=group_quantiles;
    Sample(kk).group_iqr=group_iqr;
    Sample(kk).group_lower_median=group_lower_median;
    Sample(kk).group_upper_median=group_upper_median;
    Sample(kk).group_max=group_max;
    Sample(kk).group_min=group_min;
    Sample(kk).group_sem=group_sem;
    Sample(kk).group_mean_rel=group_mean_rel;
    Sample(kk).group_std_rel=group_std_rel;
    Sample(kk).group_sem_rel=group_sem_rel;
end
%% compute statistics for each whole sample

for kk=1:3

    sample_mean = mean(Sample(kk).MpClean,1);
    sample_std  = std(Sample(kk).MpClean,1);
    sample_median = median(Sample(kk).MpClean,1);
    sample_quantiles = quantile(Sample(kk).MpClean,[.25,.75],1);
    sample_iqr = sample_quantiles(2,:)-sample_quantiles(1,:);
    sample_lower_median = sample_median - 1.57.*(sample_iqr)./sqrt(numel(Sample(kk).GpClean));
    sample_upper_median = sample_median + 1.57.*(sample_iqr)./sqrt(numel(Sample(kk).GpClean));
    
    outlierIndices=zeros(numel(Sample(kk).GpClean),5); %initialize binary array to track outlier indices
    for ii=1:5
        for jj=1:numel(Sample(kk).GpClean)
            if (Sample(kk).MpClean(jj,ii)>sample_quantiles(2,ii)+1.5*sample_iqr(ii))
                outlierIndices(jj,ii)=1;
            end
            if (Sample(kk).MpClean(jj,ii)<sample_quantiles(1,ii)-1.5*sample_iqr(ii))
                outlierIndices(jj,ii)=1;
            end
        end
    end
    
    % replace outlier values with NaN
    Mp_noOutliers_sample=Sample(kk).MpClean;
    Mp_noOutliers_sample(find(outlierIndices==1))=NaN;
    
    %calculate sample max and min without outlier values removed
    sample_max=max(Mp_noOutliers_sample,[],1);
    sample_min=min(Mp_noOutliers_sample,[],1);

    % Write key variables to sample struct
    Sample(kk).sample_mean=sample_mean;
    Sample(kk).sample_std=sample_std;
    Sample(kk).sample_median=sample_median;
    Sample(kk).sample_quantiles=sample_quantiles;
    Sample(kk).sample_iqr=sample_iqr;
    Sample(kk).sample_lower_median=sample_lower_median;
    Sample(kk).sample_upper_median=sample_upper_median;
    Sample(kk).sample_max=sample_max;
    Sample(kk).sample_min=sample_min;

end
%% Plot
close all

% Create plot one row at a time, one for each variable of interest: CLAA,
% 121_size, a lattice paramater, c lattice paramater. Each row consists of
% a map of the parameter for each sample (3x), and a box plot containing
% comparisons of the 2 clusters for each sample on a shared abscissa

fig=figure;
set(fig,'Units','inches','Position',fig_pos)

% manually label clusters. Note: this may change as k-means 
% clustering assigns cluster number 
% cats = categorical(["RH","RT","IR"]);
% cats = reordercats(cats,["RH","RT","IR"]);
cats = categorical(["RH","IR"]);
cats = reordercats(cats,["RH","IR"]);

cbarxoffset = 0;
cbaryoffset = 0;

for ii=1:4
    %loop over samples
    for kk=1:3
        pan_ind = (ii-1)*4+kk;
        xv = [0,diff(Sample(kk).x_range([1,end]))*dx];
        yv = [0,(size(Sample(kk).M,1)-1)*dy];
        ax(pan_ind) = subplot(5,17,((ii-1)*17)+(kk-1)*4+[1:4]);
        if kk==1
            imagesc(yv,xv,rot90(Sample(kk).M(:,:,ii),3));
        else
            imagesc(xv,yv,Sample(kk).M(:,:,ii))
        end
        text(0,1.15,panl{pan_ind},'Units','normalized','FontSize',FontSize_Panel,'Color','k');
        ax(kk).TickDir = 'out';
        colormap(ax(pan_ind),cmapstr{ii,kk});
        cbar(pan_ind) = colorbar;
        cbar(pan_ind).FontSize = FontSize_Tick;
        cbar(pan_ind).TickDirection = 'out';
        ap = ax(pan_ind).Position;
        if ii<5   
%             cbar(pan_ind).Position = [ap(1)+ap(3)+cbarxoffset,cbaryoffset,0.005,0.15];
            text(ax(pan_ind),1.4,1.10,cbar_str{ii},...
                'Units','normalized',...
                'FontSize',FontSize_Tick,...
                'Color','k',...
                'HorizontalAlignment','center');
        end
       
        daspect([yAsp,xAsp,1])
        ax(pan_ind).FontSize = FontSize_Tick;
        xlabel(XlabelStr,'FontSize',FontSize_Label);
        ax(pan_ind).TickDir = 'out';
        ax(pan_ind).CLim = climsSample{kk,ii};
        
    end
    linkaxes([ax((ii-1)*4+1),ax((ii-1)*4+2),ax((ii-1)*4+3)],'xy')
end

cats = categorical(["RH","IR/RT"]);
cats = reordercats(cats,["RH","IR/RT"]);

% plot Cluster maps (old style)
% ii=5;
% for kk=1:3
%     pan_ind = (ii-1)*4+kk;
%     xv = [0,diff(Sample(kk).x_range([1,end]))*dx];
%     yv = [0,(size(Sample(kk).M,1)-1)*dy];
%     ax(pan_ind) = subplot(5,17,((ii-1)*17)+(kk-1)*4+[1:4]);
%     if kk==1
%         imagesc(yv,xv,rot90(Sample(kk).M(:,:,5),3));
%     else
%         imagesc(xv,yv,Sample(kk).M(:,:,5))
%     end
%     text(0,1.15,panl{pan_ind},'Units','normalized','FontSize',FontSize_Panel,'Color','k');
%     ax(pan_ind).TickDir  = 'out';
%     colormap(ax(pan_ind),'gray');
%     cbar(pan_ind) = colorbar;
%     cbar(pan_ind).FontSize = FontSize_Tick;
%     cbar(pan_ind).TickDirection = 'out';
%         
%     ax(pan_ind).CLim = climsSample{kk,5};
%     daspect([yAsp,xAsp,1])
%     ax(pan_ind).FontSize = FontSize_Tick;
%     xlabel(XlabelStr,'FontSize',FontSize_Label);
%     
%     
%     ax(pan_ind+4)=axes('Parent',fig,...
%                'Units',ax(pan_ind).Units,'Position',ax(pan_ind).Position,...
%                'Color','none');
%     
%     % determine transparency map 
%     tmap = ones(size(Sample(kk).G));
%     % set colormap
%     % cmap = jet(3);
%     cmap = [0 0 1;0 .75 0];
%     
%     % superpose transparent cluster map
%     alpha_val = 0.15;
%     if kk==1
%         imagesc(ax(pan_ind+4),yv,xv,rot90(Sample(kk).G,3),'alphadata',alpha_val*rot90(tmap,3));
%     else    
%         imagesc(ax(pan_ind+4),xv,yv,Sample(kk).G,'alphadata',alpha_val*tmap);
%     end
%     daspect([yAsp,xAsp,1])
%     colormap(ax(pan_ind+4),cmap);
%     caxis(ax(pan_ind+4),[0.5,3.5]);
%     ax(pan_ind+4).Visible = 'off';
%     cbar(pan_ind).Visible = 'off';
% 
%     linkprop([ax(pan_ind),ax(pan_ind+4)],'Position');
%     linkaxes([ax(pan_ind),ax(pan_ind+4)]);
%     % draw and format colorbar
%     cbar(pan_ind+4) = colorbar;
%     cbar(pan_ind+4).Ticks  = 1:2:3;
%     cbar(pan_ind+4).TickLabels = {'RH','IR/RT'};
%     cbar(pan_ind+4).Direction = 'reverse';
%     cbar(pan_ind+4).TickDirection = 'out';
%     cbar(pan_ind+4).FontSize = FontSize_Tick;
%     ap = ax(pan_ind+4).Position;
%     cbar(pan_ind+4).Position = [ap(1)+ap(3)+0.005,.15,0.005,0.05];
% 
% end
% % linkaxes([ax(1),ax(2),ax(3),ax(4),ax(5),ax(6)],'xy')

% plot Cluster maps (rgb style)
ii=5;
for kk=1:3
    pan_ind = (ii-1)*4+kk;
    xv = [0,diff(Sample(kk).x_range([1,end]))*dx];
    yv = [0,(size(Sample(kk).M,1)-1)*dy];
    ax(pan_ind) = subplot(5,17,((ii-1)*17)+(kk-1)*4+[1:4]);
    
    ClusterMap_FalseColor = zeros(size(Sample(kk).M,1),size(Sample(kk).M,2),3);
    ClusterMap_FalseColor(:,:,3) = rescale((Sample(kk).G==1).*Sample(kk).M(:,:,5));
    ClusterMap_FalseColor(:,:,2) = rescale((Sample(kk).G==2).*Sample(kk).M(:,:,5))+.9*rescale((Sample(kk).G==1).*Sample(kk).M(:,:,5));    
        
    if kk==1
        imagesc(rot90(ClusterMap_FalseColor,3))
    else
        imagesc(ClusterMap_FalseColor)
    end
    text(0,1.15,panl{pan_ind},'Units','normalized','FontSize',FontSize_Panel,'Color','k');
    ax(pan_ind).TickDir  = 'out';
    cbar(pan_ind) = colorbar;
    cbar(pan_ind).FontSize = FontSize_Tick;
    cbar(pan_ind).TickDirection = 'out';       
    ax(pan_ind).CLim = climsSample{kk,1};
    daspect([yAsp,xAsp,1])
    ax(pan_ind).FontSize = FontSize_Tick;
    xlabel(XlabelStr,'FontSize',FontSize_Label);
    ax(pan_ind+4)=axes('Parent',fig,...
               'Units',ax(pan_ind).Units,'Position',ax(pan_ind).Position,...
               'Color','none');
    daspect([yAsp,xAsp,1])    

    cmap = [0 1 1;0 1 0];

    colormap(ax(pan_ind+4),cmap);
    caxis(ax(pan_ind+4),[0.5,3.5]);
    ax(pan_ind+4).Visible = 'off';
    cbar(pan_ind).Visible = 'off';

    linkprop([ax(pan_ind),ax(pan_ind+4)],'Position');
    linkaxes([ax(pan_ind),ax(pan_ind+4)]);
    % draw and format colorbar
    cbar(pan_ind+4) = colorbar;
    cbar(pan_ind+4).Ticks  = 1:2:3;
    cbar(pan_ind+4).TickLabels = {'RH','IR/RT'};
    cbar(pan_ind+4).Direction = 'reverse';
    cbar(pan_ind+4).TickDirection = 'out';
    cbar(pan_ind+4).FontSize = FontSize_Tick;
    ap = ax(pan_ind+4).Position;
    cbar(pan_ind+4).Position = [ap(1)+ap(3)+0.005,.15,0.005,0.05];

end

linkaxes([ax((ii-1)*4+1),ax((ii-1)*4+2),ax((ii-1)*4+3)],'xy')


%% boxplots for panels D, H, L, P, T

% Assemble sample datasets into singular matrix with 2 grouping variables:
% sample, and RH vs RT/IR
cats = categorical(["RH","IR/RT"]);
cats = reordercats(cats,["RH","IR/RT"]);
samples = categorical(["S1","S2","S3"]);

MpCleanAll = [Sample(1).MpClean;Sample(2).MpClean;Sample(3).MpClean];
GpCleanAll = [Sample(1).GpClean,ones(numel(Sample(1).GpClean),1);Sample(2).GpClean,2*ones(numel(Sample(2).GpClean),1 );Sample(3).GpClean,3*ones(numel(Sample(3).GpClean),1)];
GpCleanAllCats=[samples(GpCleanAll(:,2))',cats(GpCleanAll(:,1))'];

for ii=1:4
    pan_ind=ii*4;
    ax(pan_ind) = subplot(5,17,((ii-1)*17)+[14:17]);
    boxplot(MpCleanAll(:,ii),GpCleanAllCats,...
            'Notch','on',...
            'BoxStyle','outline',...
            'DataLim',dlim_abs(ii,:),...
            'ExtremeMode','clip',...
            'Symbol','.r',...
            'Widths',.5,...
            'FactorGap',[10,2],...
            'OutlierSize',2,...
            'LabelVerbosity','minor');
            hold on;
    ylim(climsGlobal{ii});
%     set(gca,'XTickLabel',{' '})
%     % fill boxes using colors from group map
%     % note that the order in h is the reverse of the axis for some reason
    h = findobj(ax(pan_ind),'Tag','Box');
    for kk=1:length(h)
         patch(get(h(kk),'XData'),get(h(kk),'YData'),cmap(mod(kk,2)+1,:),'EdgeColor',cmap(mod(kk,2)+1,:));
    end

%     % annotate with group mean
%     plot(group_mean(:,ii),'r_')
    % format plot
    ax(pan_ind).YMinorTick = 'on';
    ax(pan_ind).TickDir = 'out';
    ax(pan_ind).FontSize = FontSize_Tick;
% 
%     % indicate mean and standard deviation across all clusters
%     yline(mean_val(ii),'k-','LineWidth',1);
%     yline(mean_val(ii)+std_val(ii),'k:','LineWidth',1);
%     yline(mean_val(ii)-std_val(ii),'k:','LineWidth',1);
    xlabel('Group','FontSize',FontSize_Label)
    ylabel(join(cbar_str{ii}),'FontSize',FontSize_Label);
    text(ax(pan_ind),-0.15,1.1,panl{pan_ind},'FontSize',FontSize_Panel,'Color','k','Units','normalized')
%     yyaxis right
%     ax(ii+5).YAxis(2).Limits=100*(ax(ii+5).YAxis(1).Limits/mean_val(ii)-1);
%     ax(ii+5).YAxis(2).MinorTick = 'on';
%     ylabel(join(cbar_str2{ii}),'FontSize',FontSize_Label)
    hold off
end

%% Save figure
if save_figure
    exportgraphics(gcf,[pn_out,fn_out{1}],'BackgroundColor','none','ContentType','vector');
    disp(['Saved as: ',pn_out,fn_out{1}])
else
    disp('Figure not saved!');
end
    


if save_tables
    writetable(T1,[pn_out,fn_out{3}]);
    writetable(T2,[pn_out,fn_out{2}]);
    disp(['Saved as: ',pn_out,fn_out{2}])
else
    disp('Figure not saved!');
end

%% Scatter Plots of Raw Data
close all
for kk=1:3
    figure
    DeltaMglist = reshape(Sample(kk).MpClean(:,3),[numel(Sample(kk).MpClean(:,3)),1]);
    DeltaCO3 = reshape(Sample(kk).MpClean(:,4),[numel(Sample(kk).MpClean(:,4)),1]);
    
    Sample(kk).CovAC = cov(DeltaMglist,DeltaCO3);
    error_ellipse(Sample(kk).CovAC,Sample(kk).sample_mean(3:4))
    hold on
    scatter(DeltaMglist,DeltaCO3);
    xlabel('\Deltax_{Mg} [at%]')
    ylabel('\Deltax_{CO_3} [at%]')
    title({'Compositional Difference Scatter: ',Sample(kk).name})
end