% Use this script with the output of Step1_ProcessEnamelDiffractionPatterns
% to perform the k-means clustering analysis presented in Fig 3.

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
save_map = true;
save_stats = true;

%% Read in data
CLAAMatrix=readmatrix('CLAA_Map.csv');
AMatrix=readmatrix('A_Parameter_Map.csv');
CMatrix=readmatrix('C_Parameter_Map.csv');
Size121Matrix=readmatrix('121_Scherrer_Size_Map.csv');

%% Assemble data for cluster analysis
%Define dimensions of map in pixels (different for each scan)
xPts=52;
yPts=15;

%number of patterns
numPatt=xPts*yPts;
pattRange=1:numPatt;
pattMatrix=reshape(pattRange,[xPts,yPts])';
% pattSampleMatrix=pattMatrix(:,:); %Cropped area based on examining results
% pattSample=reshape(pattSampleMatrix',[1,numel(pattSampleMatrix)]); %Create list of patterns determined to be sample
% pattPlatinum=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,38,39,40,41,42,43,44,45,46,47,48,49,50,51,75,76,77,78,79,80,81,82,83,84,85,112,113,114,115,149,150,151,186,187,778,815,852,889,890]; %the number index of each pattern identified to have significant platinum diffraction
% indexMatrix=reshape(pattRange,[yPts,xPts])'; %create a matrix oriented such that indexing by pattern # will yield the corresponding matrix element index
% PlatIndices=indexMatrix(pattPlatinum); %matrix indices that can be passed to yield an output from patterns associated with high Pt patterns

%set Platinum points to zero
% CLAAMatrix(PlatIndices)=0;
% AMatrix(PlatIndices)=0;
% CMatrix(PlatIndices)=0;
% Size121Matrix(PlatIndices)=0;

D(:,:,1) = CLAAMatrix;
D(:,:,2) = AMatrix;
D(:,:,3) = CMatrix;
D(:,:,4) = Size121Matrix;

varname = {'CLAA','a','c','s_{121}'};

I = D(:,4:52,:); % crop patterns from left side
Ip = squeeze(reshape(I,[],1,4)); % reshape to vector

samplePixels = find(Ip(:,1)>0); % find nonzero pixels where good sample is
IpClean=Ip(samplePixels,:);
zIpClean=zscore(IpClean);

%reassemble the full list with NaNs, but with the zscore calculated
%without zeros.
zIp=Ip;
jj=1;
for ii=1:length(Ip)
    if Ip(ii,1)>0
        zIp(ii,:)=zIpClean(jj,:);
        jj=jj+1;
    else
        zIp(ii,:)=[NaN,NaN,NaN,NaN];
    end    
end

%% k-means clustering all parameters
k = 2;
dist = 'correlation';
[GpClean,centr3,sumdist3] = kmeans(zIpClean(:,[1,2,3,4]),k,'Distance',dist,'Display','final','Replicates',10);

% reassemble with NaNs
Gp=NaN(length(Ip),1);
jj=1;
for ii=1:length(Gp)
    if Ip(ii,1)>0
        Gp(ii)=GpClean(jj);
        jj=jj+1;
    else
        Gp(ii)=NaN;
    end    
end

% reshape to matrix
G = reshape(Gp,15,[]);

figure(1);
[silh3,h] = silhouette(zIp(:,[1,2,3,4]),Gp,dist);
xlabel('Silhouette Value')
ylabel('Cluster')
title('Cluster Silhouettes as Quality Metric')

% determine averages using cluster ID as grouping variable
means_by_group_kmeans = groupsummary(Ip,Gp,'mean');
std_by_group_kmeans   = groupsummary(Ip,Gp,'std');

figure(2);
clear ax
for ii = 1:4
    ax{ii} = subplot(1,4,ii);
    bar(means_by_group_kmeans(:,ii))
    ax{ii}.TickDir = 'out';
    ylabel(varname{ii})
    xlabel('cluster ID')
    hold on
    er = errorbar(means_by_group_kmeans(:,ii),std_by_group_kmeans(:,ii));
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    hold off
end
sgtitle('Mean Crystallographic Parameters by Cluster') 

ax{1}.YLim = [0,8];
ax{2}.YLim = [9.442,9.456];
ax{3}.YLim = [6.85,6.91];
ax{4}.YLim = [350,550];

%% k-means clustering no crystallite size
k = 2;
dist = 'correlation';
[Gp,centr3,sumdist3] = kmeans(zIp(:,[1,2,3]),k,'Distance',dist,'Display','final','Replicates',10);

% reshape to matrix
G = reshape(Gp,15,[]);

figure(1);
[silh3,h] = silhouette(zIp(:,[1,2,3]),Gp,dist);
xlabel('Silhouette Value')
ylabel('Cluster')
title('Cluster Silhouettes as Quality Metric')

% determine averages using cluster ID as grouping variable
means_by_group_kmeans = groupsummary(Ip,Gp,'mean');
std_by_group_kmeans   = groupsummary(Ip,Gp,'std');

figure(2);
clear ax
for ii = 1:3
    ax{ii} = subplot(1,3,ii);
    bar(means_by_group_kmeans(:,ii))
    ax{ii}.TickDir = 'out';
    ylabel(varname{ii})
    xlabel('cluster ID')
    hold on
    er = errorbar(means_by_group_kmeans(:,ii),std_by_group_kmeans(:,ii));
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    hold off
end
sgtitle('Mean Crystallographic Parameters by Cluster') 

ax{1}.YLim = [0,8];
ax{2}.YLim = [9.442,9.456];
ax{3}.YLim = [6.85,6.91];

%% overlay clusters
figure(3);
ax13 = axes;
imagesc(I(:,:,1));
ax13.TickDir = 'out';
axis image
colormap(ax13,'gray');
ax23 = axes;


tmap = zeros(size(G));
for jj=1:k
    tmap = tmap | G==jj;
end

imagesc(ax23,G,'alphadata',0.25*tmap);
axis image
colormap(ax23,jet(max(Gp)));
caxis(ax23,[min(nonzeros(G)), max(nonzeros(G))]);
ax23.Visible = 'off';
linkprop([ax13 ax23],'Position');
cbar = colorbar;
cbar.Ticks = sort(1:k);
cbar.TickDirection = 'out';

%plot just cluster assignments
figure
imagesc(G)
daspect([1,1,1])
title('Cluster Map')

%% 1-way ANOVAs
figure(1)
[p,tbl,stats] = anova1(Ip(:,3),Gp,'off')
[c,m,h,gnames] = multcompare(stats)
[gnames num2cell(m)]

clc
TS3 = table('Size',[0,7],...
          'VariableTypes',{'categorical','categorical','categorical','double','double','double','double'},...
          'VariableNames',{'Varname','Group1','Group2','LowerCI','Difference','UpperCI','p'});
      
for ii = 1:size(Ip,2)
    [~,~,stats] = anova1(Ip(:,ii),Gp,'off');
    [c,m,~,gnames] = multcompare(stats,'Display','off');
    Temp = array2table(c,'VariableNames',{'Group1','Group2','LowerCI','Difference','UpperCI','p'});
    Temp.Varname = categorical(cellstr(repmat(varname{ii},height(Temp),1)));
    Temp.Group1 = categorical(gnames(Temp.Group1));
    Temp.Group2 = categorical(gnames(Temp.Group2));
    TS3 = [TS3;Temp];
end
disp(TS3)



%% write results
if save_map
    writematrix(G,'Cluster_Map.csv');
else
    disp('Cluster Map Not Saved!')
end

if save_stats
    writetable(TS3,'MultipleComparison_TS3.csv');
else
    disp('Stats not Saved!')
end
