% Use this script in conjunction with the .csv data files to generate rough
% versions of the figures and tables for 
% "Mesoscale Structural Gradients in Human Tooth Enamel", by R. Free, 
% K. DeRocher, V. Cooley, R. Xu, S.R. Stock, and D. Joester.
% This script also includes the units and axes information for each plot.

% Author: Robert Free and Derk Joester

%% Clean Up
clear variables 
close all
clc

%% Flags
save_figs = false;

%% Load Data
mfile_name          = 'Fig5_and_S7_plot.m';

% paths & filenames
pn_data = './figure source data/';
fn_inp     = {'ContourMap_from_Boyde.csv',...
              'Human_Enamel_TPP_Boyde_1.png'};

pn_out = './';
          
% Import data
M = readmatrix([pn_data,fn_inp{1}],'NumHeaderLines',3);
xv = M(1,2:end);
yv = M(2:end,1);
T  = M(2:end,2:end);
I = imread([pn_data,fn_inp{2}]);

%% define scale and offset for original image and digitized contours
sx = 62.2; % pixels/µm
sy = 61; % pixels/µm
offx = 11/sx;

% Define plot options
figS7_width     = 7;
figS7_asp_ratio = 28/18;
figS7_zoom      = 4;
figS7_pos       = [3,8,1+figS7_width*figS7_zoom,figS7_width*figS7_zoom/figS7_asp_ratio]; %[in]

FontSize_Panel = 10*figS7_zoom;
FontSize_Label = 7*figS7_zoom;
FontSize_Tick  = 6*figS7_zoom;

%% calculate contour matrix
CM = contourc(xv,yv,T,unique(T(~isnan(T))));

%% convert contour matrix to point cloud
m(1)=1; 
n=1;  
try
    while n<length(CM)
        n=n+1;
        m(n) = m(n-1)+CM(2,m(n-1))+1;         
    end
end
XYZ = [NaN,NaN,NaN]; 
for nn = 1:n-2
     x{nn} = CM(1,m(nn)+1:m(nn+1)-1); 
     y{nn} = CM(2,m(nn)+1:m(nn+1)-1);    
     z(nn) = CM(1,m(nn));
     XYZ = [XYZ;[x{nn}',y{nn}',repmat(z(nn),length(x{nn}),1)]];
end
XYZ = XYZ(2:end,:);

% fit plane to point cloud
b         = [XYZ(:,1:2),ones(length(XYZ),1)]\XYZ(:,3);
% determine normal vector of fitted plane
n = [b(1);b(2);-1];
n_hat = n/norm(n);

clear m x y z n
%% plot
[INx, INy, ~] = size(I);

% define plot box for 3D plot
x_min = 0;
x_max = 20;
y_min = 0;
y_max = 20;
z_min = 0;
z_max = 25;
Lims = [[x_min,x_max];[y_min,y_max];[z_min,z_max]];
% define edges of plot box
dh = [[1;0;0],[0;1;0],[0;0;1]];
p  = [[x_min;y_min;z_min],[x_min;y_max;z_min],[x_min;y_min;z_max],[x_min;y_max;z_max],...
      [x_min;y_min;z_min],[x_max;y_min;z_min],[x_min;y_min;z_max],[x_max;y_min;z_max],...
      [x_min;y_min;z_min],[x_max;y_min;z_min],[x_min;y_max;z_min],[x_max;y_max;z_min]];
  
  
figS7 = figure('Units','inches','Position',figS7_pos);
ax1 = subplot(231);
imagesc([0,(INx-1)/sx], [0,(INy-1)/sy], sum(I,3)); hold on;
colormap(ax1,'gray');
axis image
ax1.TickDir = 'out';
ax1.FontSize = FontSize_Tick;
ax1.XMinorTick = 'on';
ax1.YMinorTick = 'on';
xlim([0,17])
xlabel('x [µm]','FontSize',FontSize_Label);
ylabel('y [µm]','FontSize',FontSize_Label);
text(ax1,-3,0,'A', 'FontSize',FontSize_Panel,'VerticalAlignment','middle');
ax1.XAxis.LineWidth = 2;
ax1.YAxis.LineWidth = 2;

ax2 = subplot(232);
imagesc(xv,yv,T);
colormap(ax2,'parula')
axis image
ax2.TickDir = 'out';
ax2.XMinorTick = 'on';
ax2.YMinorTick = 'on';
ax2.FontSize = FontSize_Tick;
cbar2 = colorbar;
cbar2.TickDirection = 'out';
cbar2.Label.String = 'height [µm]';
xlabel('x [µm]');
ylabel('y [µm]');
text(ax2,-3,0,'B', 'FontSize',FontSize_Panel,'VerticalAlignment','middle');
ax2.XAxis.LineWidth = 2;
ax2.YAxis.LineWidth = 2;

ax3 = subplot(233);
plot(XYZ(:,1)+offx,XYZ(:,2),'r.')
axis image
ax3.TickDir = 'out';
ax3.FontSize = FontSize_Tick;
ax3.YDir = 'reverse';
ax3.XMinorTick = 'on';
ax3.YMinorTick = 'on';
xlim([0,17])
xlabel('x [µm]','FontSize',FontSize_Label);
ylabel('y [µm]','FontSize',FontSize_Label);
text(ax3,-3,0,'C', 'FontSize',FontSize_Panel,'VerticalAlignment','middle');
ax3.XAxis.LineWidth = 2;
ax3.YAxis.LineWidth = 2;

ax4 = subplot(234);
imagesc([0,(INx-1)/sx], [0,(INy-1)/sy], sum(I,3)); hold on;
colormap(ax4,'gray');
plot(XYZ(:,1)+offx,XYZ(:,2),'r.')
axis image
ax4.TickDir = 'out';
ax4.XMinorTick = 'on';
ax4.YMinorTick = 'on';
ax4.FontSize = FontSize_Tick;
xlim([0,17])
xlabel('x [µm]','FontSize',FontSize_Label);
ylabel('y [µm]','FontSize',FontSize_Label);
text(ax4,-3,0,'D', 'FontSize',FontSize_Panel,'VerticalAlignment','middle');
ax4.XAxis.LineWidth = 2;
ax4.YAxis.LineWidth = 2;

% plot point cloud and fitted plane
ax5 = subplot(235);
plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'.');hold on;
grid on;
ctr  = 0;
Vert = [];
    for jj=1:3
        for kk = 1:4
            ctr = ctr+1;
            y = intersect_plane_line([0;0;b(3)],n_hat,dh(:,jj),p(:,ctr),1e-9);
            if (~isempty(y) | ~isnan(y)) & (y(jj) >= Lims(jj,1)) & (y(jj) <= Lims(jj,2))
                Vert = [Vert,y];
                Face = [1,2,4,3];
            end
        end
    end
patch(ax5,'Vertices',Vert','Faces',Face,'FaceColor','r','FaceAlpha',0.3,'EdgeColor','none');
xlabel('x [µm]','FontSize',FontSize_Label);
ylabel('y [µm]','FontSize',FontSize_Label);
zlabel('z [µm]','FontSize',FontSize_Label);
ax5.TickDir = 'out';
ax5.FontSize = FontSize_Tick;
ax5.XMinorTick = 'on';
ax5.YMinorTick = 'on';
ax5.YMinorTick = 'on';
text(ax5,-1,24,25,'E', 'FontSize',FontSize_Panel,'VerticalAlignment','middle');
ax5.XAxis.LineWidth = 2;
ax5.YAxis.LineWidth = 2;
ax5.ZAxis.LineWidth = 2;
% prepare panel 
ax6 = subplot(236);
clear y Vert Face

%% create surface model from point cloud
% define scattered Interpolant
uXYZ = unique(XYZ,'rows');
F = scatteredInterpolant(uXYZ(:,1),uXYZ(:,2),uXYZ(:,3),'natural','none');

% define grid vectors
x2v =0:0.1:18; 
y2v = 0:0.1:20;
% define grid mesh
[X2g,Y2g] = meshgrid(x2v,y2v);
[Nx,Ny] = size(X2g);

% interpolate over mesh
Z2g = F(X2g,Y2g);

% angle between normal vector of fitted plane and [001]
qt = vrrotvec(n_hat,[0;0;1]);
% rotation matrix
R1 = vrrotvec2mat(qt);

% translation
XYZ2 = [X2g(:),Y2g(:),Z2g(:)];
Tv = mean(XYZ2,'omitnan');

% rotate point cloud to align normal vector of fitted plane with [001],
% then rotate to align edges with [100] and [010]
R2 = vrrotvec2mat([0,0,1,pi/9]);

XYZ2r = (R2*R1*(XYZ2-Tv)')'; % 36381-by-3
XYZ3r = reshape(XYZ2r,Nx,Ny,3); % Nx-by-Ny-by-3

% determine surface normals
[nx,ny,nz] = surfnorm(XYZ3r(:,:,1)',XYZ3r(:,:,2)',imgaussfilt(XYZ3r(:,:,3),2)','EdgeColor','none');

% plot rotated point cloud with smoothing
fig5A = figure;
ax7 = axes;
surf(ax7,XYZ3r(:,:,1),XYZ3r(:,:,2),imgaussfilt(XYZ3r(:,:,3)-min(XYZ3r(:,:,3)+1,[],'all'),2),...
    'EdgeColor','none',...
    'FaceColor','interp',...
    'FaceLighting','gouraud', ...
    'AmbientStrength',0.5,...
    'DiffuseStrength',0.7, ...
    'SpecularStrength',0.5,...
    "LineWidth",2);

hold on;
daspect([1 1 1])
xlim([-10,10]);
ylim([-10,10]);
zlim([-1,5]);

view([-148 22])
camlight
xlabel('x [µm]')
ylabel('y [µm]')
zlabel('z [µm]')
cb1 = colorbar;
cb1.Position = [0.92 0.500 0.02 0.2000];
cb1.TickDirection = 'out';
cb1.Label.String = 'height';
colormap turbo
hold off

if save_figs
    saveas(fig5A,[pn_out,'Fig5A.eps'],'epsc');
end

%% Figure 5B: Profiles of Tomes's Processes

% define initial unit normal vector for all planes (xy plane rotated by 60˚ cw about [100])
P = [[0,0,0];...
     [10,0,0];
     [0,0.5,sqrt(3)/2]]';

n2    = cross(P(:,2)-P(:,1),P(:,3)-P(:,1));
n2hat = -n2/norm(n2);

% define rotation of unit normal vector about the "vertical" direction in
% the initial plane
raxis = cross(cross(n2hat,[0;0;1]),n2hat);
ang   = [-pi/36,-pi/45,-pi/72,-pi/48,0]; % angles picked by trial and error

% manually define one point on each plane
Q = [[-5.546,-4.394,2.744];...
     [2.273,-2.8,1.559];...
     [-5.298,1.145,1.915]
     [0.4679,4.538,2.209];...
     [-8.251,7.803,3.097]]';

% some preliminary calcs for intersections of planes with surface
DT = delaunay(X2g,Y2g);
tol = 1e-9;
X3rg = XYZ3r(:,:,1);
Y3rg = XYZ3r(:,:,2);
Z3rsg = imgaussfilt(XYZ3r(:,:,3),2);

% for coloring, determine cosine of angle between a manually selected
% reference direction and the surface normal at the surface vertex
ref_az = 69.11; % [˚]
ref_el = -28;   % [˚]

[refx,refy,refz] = sph2cart(ref_az,ref_el,1);
C_cos = rescale(refx*nx + refy*ny + refz*nz)';

% define plot box
x_min = -10;
x_max = 10;
y_min = -10;
y_max = 10;
z_min = -3;
z_max = 5;
Lims = [[x_min,x_max];[y_min,y_max];[z_min,z_max]];

surf(ax6,X3rg,Y3rg,Z3rsg,repmat(C_cos,1,1,3),'EdgeColor','none','FaceColor','interp');
hold(ax6,'on');
text(ax6,16,10,5,'F', 'FontSize',FontSize_Panel,'VerticalAlignment','middle');
daspect(ax6,[1 1 1])
lighting(ax6,'none');
xlim(ax6,[x_min,x_max]);
ylim(ax6,[y_min,y_max]);
zlim(ax6,[z_min,z_max]);

view(ax6,[-114.74 31.76]);
xlabel(ax6,'x [µm]', 'FontSize',FontSize_Label)
ylabel(ax6,'y [µm]', 'FontSize',FontSize_Label)
zlabel(ax6,'z [µm]', 'FontSize',FontSize_Label)
caxis(ax6,[0.4,.6])
ax6.TickDir = 'out';
ax6.XMinorTick = 'on';
ax6.YMinorTick = 'on';
ax6.YMinorTick = 'on';
ax6.FontSize = FontSize_Tick;
ax6.XAxis.LineWidth = 2;
ax6.YAxis.LineWidth = 2;
ax6.ZAxis.LineWidth = 2;
%% plot profiles superimposed on surface (ax1) and in plane of intersection (ax2)
fig5B = figure;
ax8 = axes;
ax8.YMinorTick = 'on';
ax8.XMinorTick = 'on';
ax8.TickDir = 'out';
daspect([1,1,1])
grid on
xlabel('distance [µm]')
ylabel('distance [µm]')
xlim([-15,5]);
ylim([-1,26]);
% plot parameters for line profiles in lane of intersection
offset = -5;
colm = parula(6);

% define fit boundaries
fit_lb = -7;
fit_ub = -0.5;

% define plot boundaries for fitted curve
plot_lb = -8;
plot_ub = 1;

% define edges of plot box
dh = [[1;0;0],[0;1;0],[0;0;1]];
p  = [[x_min;y_min;z_min],[x_min;y_max;z_min],[x_min;y_min;z_max],[x_min;y_max;z_max],...
      [x_min;y_min;z_min],[x_max;y_min;z_min],[x_min;y_min;z_max],[x_max;y_min;z_max],...
      [x_min;y_min;z_min],[x_max;y_min;z_min],[x_min;y_max;z_min],[x_max;y_max;z_min]];

for ii=1:5
    % rotate plane unit normal vector
    n2hr = vrrotvec2mat([raxis;ang(ii)])*n2hat;
    
    % determine intersection of plane defined by n2hr and a point Q with
    % smoothed surface
    [XYZv,conn,XYp] = intersect_plane_surf(Q(:,ii),n2hr,X3rg,Y3rg,Z3rsg,DT,tol);
    
%     % determine vertices for patch object to represent rotated plane
%     % direction of intersection of rotated plane with xy plane
%     pv1 = cross(n2hr,[0;0;1]);
%     % direction normal to both n2hr and pv1
%     pv2 = cross(pv1,n2hr);
    clear y
    ctr  = 0;
    Vert = [];
    for jj=1:3
        for kk = 1:4
            ctr = ctr+1;
            y = intersect_plane_line(Q(:,ii),n2hr,dh(:,jj),p(:,ctr),1e-9);
            if (~isempty(y) | ~isnan(y)) & (y(jj) >= Lims(jj,1)) & (y(jj) <= Lims(jj,2))
                Vert = [Vert,y];
                Face = [1,2,4,3];
            end
        end
    end
    patch(ax6,'Vertices',Vert','Faces',Face,'FaceColor',colm(ii,:),'FaceAlpha',0.3,'EdgeColor','none');
    % each intersection consists of a sequence of lines that result from
    % the intersection of a triangle with the plane, so plot them in
    % sequence
    for kk = 1:length(conn)
        % plot on top of surface
        plot3(ax6,XYZv(conn(kk,:),1),XYZv(conn(kk,:),2),XYZv(conn(kk,:),3),'-','Color',colm(ii,:),'LineWidth',2); hold on;
    end
    % in principle, an intersection could consist of more than one
    % connected component, but here we selected the planes such that that
    % doesn't happen.
    for jj = 1:1
        % determine indices for fitting ...
        ind = (XYp{jj}(1,:)>=fit_lb) & (XYp{jj}(1,:)<=fit_ub);
        % and plotting of fit
        ind2 = (XYp{jj}(1,:)>plot_lb) & (XYp{jj}(1,:)<plot_ub);
        % fit linear function to "ramp" of TP
        f{ii} = fit(XYp{jj}(1,ind)',XYp{jj}(2,ind)'+ii*offset+30,'poly1');
        slope(ii) = f{ii}.p1;
        
        % plot
        plot(ax8,XYp{jj}(1,:),XYp{jj}(2,:)+ii*offset+30,'.','Color',colm(ii,:),'MarkerSize',2); hold on;
        plot(ax8,XYp{jj}(1,ind2),f{ii}(XYp{jj}(1,ind2)),'r:');
        
    end
end
xline(ax8,fit_lb,'k:');
xline(ax8,fit_ub,'k:');

mean_slope_angle = mean(atand(slope));
std_slope_angle = std(atand(slope));

text(ax8,-5,25,{['$\bar{\alpha} = ',num2str(mean_slope_angle,'%2.1f'),'^\circ $'],...
                ['$\sigma_\alpha = ',num2str(std_slope_angle,'%2.1f'),'^\circ $']},...
                'HorizontalAlignment','center','interpreter','latex');
view(ax6,-114.4,21.4);
if save_figs
    saveas(figS7,[pn_out,'FigS7.pdf'],'pdf');
    saveas(fig5B,[pn_out,'Fig5B.eps'],'epsc');
end

%% Figure 5D,E: TPs and pits shaded by "age" i.e. distance from section plane to growth front

XYZ3rs(:,:,1:2) =  XYZ3r(:,:,1:2);
XYZ3rs(:,:,3) = imgaussfilt(XYZ3r(:,:,3)-min(XYZ3r(:,:,3)+1,[],'all'),2);
XYZ2rs = reshape(XYZ3rs,[],3);

% normal vector of section plane
[nn1,nn2,nn3]=sph2cart((-99.8+90)/180*pi,-19.73/180*pi,1);
n_hat = [nn1;nn2;nn3];
q = [-11,0,0];

delta = (q - XYZ2rs)*n_hat;
Cdist = reshape(-delta,Nx,Ny);

figure;
surf(XYZ3r(:,:,1),XYZ3r(:,:,2),imgaussfilt(XYZ3r(:,:,3)-min(XYZ3r(:,:,3)+1,[],'all'),2),Cdist,...
    'EdgeColor','none',...
    'FaceColor','interp');
hold on;
daspect([1 1 1])
xlabel('x [µm]')
ylabel('y [µm]')
zlabel('z [µm]')

view([-99.8 19.73]) % view along rod direction

cb1 = colorbar;
cb1.Position = [0.92 0.500 0.02 0.2000];
cb1.TickDirection = 'out';
cb1.Label.String = 'distance [µm] (age)';
colormap(turbo(52));
caxis([-1.5,24.5]);

% scaleball with 2 µm diameter
[Xs,Ys,Zs] = sphere(36);
surf(Xs-9,Ys-3.5,Zs,...
    'EdgeColor','none',...
    'FaceColor',[0.7,0.7,0.7]);
axis off
hold off

if save_figs
    saveas(gcf,[pn_out,'Fig5E'],'epsc');
end
hold off

% same colorscheme, but show pits instead

figure;
surf(XYZ3r(:,:,1),XYZ3r(:,:,2),imgaussfilt(XYZ3r(:,:,3)-min(XYZ3r(:,:,3)+1,[],'all'),2),Cdist,...
    'EdgeColor','none',...
    'FaceColor','interp');
hold on;
daspect([1 1 1])
xlabel('x [µm]')
ylabel('y [µm]')
zlabel('z [µm]')

view([80 -64]) % view of pits

cb1 = colorbar;
cb1.Position = [0.92 0.500 0.02 0.2000];
cb1.TickDirection = 'out';
cb1.Label.String = 'distance [µm] (age)';
colormap(turbo(52));
caxis([-1.5,24.5]);

% scaleball with 2 µm diameter
[Xs,Ys,Zs] = sphere(36);
surf(Xs-10.5,Ys+3,Zs,...
    'EdgeColor','none',...
    'FaceColor',[0.7,0.7,0.7]);
axis off

lighting gouraud
camlight

hold off
if save_figs
    saveas(gcf,[pn_out,'Fig5D.eps'],'epsc');
end

%%
function y = intersect_plane_line(q,nh,dh,p,tol)
    if abs(nh(:)'*dh(:))<=tol
        if abs((p(:)-q(:))'*nh(:)) <= tol
            %exception: line is in plane
            y = [Nan];
        else
            % exception: line is parallel to plane but p does not lie on
            % plane
            y = [];
        end
    else
        x = (q(:)-p(:))'*nh(:)/(dh(:)'*nh(:));
        y = x*dh+p;
    end
end
function [XYZv,conn,XYp] = intersect_plane_surf(q,nh,Xg,Yg,Zg,DT,tol)
% intersect_plane_surf determines the intersection of a plane and a
% triangulated surface
% the plane is define by point q and unit normal vector nh
% the surface is defined by meshgrids Zg, Yg, and Zg = F(Xg,Yg)
% DT is the 2D Delaunay triangulation of (Xg,Yg)
% tol is a threshhold to identify close-by points in the intersection

    % calculate signed distance of all vertices from plane
    XYZv = [Xg(:),Yg(:),Zg(:)];
    d = (XYZv-q')*nh;
    
    % create array of signed distances for vertices of all triangles
    s = sign(d(DT));
    
    % determine the number of triangle vertices that fall onto the plane
    % todo: implemement tolerance criterion
    n = zeros(length(s),1);
    for ii = 1:length(s)
        n(ii) = 3-nnz(s(ii,:));
    end
    
    % initialize connectivity matrix
    % each row describes contains two indices of vertices in XYZv
    % the intersection is a line between the indexed vertices
    conn = [];
    
    % initialize counter for vertices added later
    ctr = length(XYZv)+1;
    
    % if all three vertices are on the plane, add 3 edges to connectivity
    % matrix
    ind3 = find(n==3);
    for ii = 1:length(ind3)
        V = DT(ind3(ii),s(ind3(ii),:)==0);
        conn = [conn;V(1:2);V([1,3]);V(2:3)];
    end
    
    % if exactly two vertices are on the plane, add 1 edge to connectivity matrix
    ind2 = find(n==2);
    for ii = 1:length(ind2)
        conn = [conn;DT(ind2(ii),s(ind2(ii),:)==0)];
    end
    
    % if one vertex is on the plane
    ind1 = find(n==1);
    for ii = 1:length(ind1)
        % and if the other two vertices are on opposite sides, find intercept of plane
        % with crossing edge
        if sum(nonzeros(s(ind1(ii),:)))==0
            % two vertices, q and r, that are on opposite side of the plane
            qr = XYZv(DT(ind1(ii),s(ind1(ii),:)~=0),:);
            % the vertex p that is on the plane
            p  = XYZv(DT(ind1(ii),s(ind1(ii),:)==0),:);
            % unit normal to the plane through the vertices
            nth = cross(p-qr(1,:),p-qr(2,:));
            nth = nth/norm(nth);
            % the unit direction vector of the line of intersection between
            % the section plane and the triangle plane
            dh1 = cross(nth,nh);
            % the unit direction vector of the triangle edge that crosses
            % the plane
            dh2 = diff(qr)/norm(diff(qr));
            % solve for intercept between the line through p with direction
            % vector dh1 and the the qr edge
            dv = [dh1;dh2]'\(qr(1,:)-p)';
            % add intercept to vertex array and increment counter
            XYZv(ctr,:) = dv(1)*dh1+p;
            ctr = ctr+1;
            % update connectivity matrix
            conn = [conn;[DT(ind1(ii),s(ind1(ii),:)==0),ctr]];           
        end
        % if the other vertices are on the same side, ignore
    end
    
    % find all remaining triangles (w/o vertices on the
    % plane) that have vertices on opposite sides of the plane
    % note that there are the following possibilities for a row in s:
    % +/-[1,1,1]; all permutations of +/-[-1,1,1]
    % if we take the absolute of the row sum, we either get 3 or 1
    % only the latter have vertices on opposite sides
    indt = find(n == 0 & abs(sum(s,2))==1);
    
    % iterate over all triangles with edges that cross the plane
    for ii = 1:length(indt)
        % grab vertices
        V = XYZv(DT(indt(ii),:),:);
        % get row vector from s for triangle
        sv = s(indt(ii),:);
        % make sure that we have a permutation of [1,-1,-1] 
        if sum(sv)==1
            sv = -1*sv;
        end
        % now vertex p is on ones side
        p  = V(sv==1,:);
        % vertices q and r on the other
        qr = V(sv==-1,:);
        
        % for the two edges that cross the plane, find the intercept of the
        % edge with the plane
        for jj=1:2
            % unit vector of direction of crossing edge
            dh = (-p+qr(jj,:))/norm(-p+qr(jj,:));
            % solve for intersection between edge and plane
            d  = ((q'-p)*nh)/(dh*nh);
            % determine new vertex of line segment
            XYZv(ctr,:) = d * dh + p;
            ctr = ctr+1;
        end
        % update connectivity matrix
        conn = [conn;[ctr-2,ctr-1]];
    end
    
    % identify unique vertices and their number
    uverts = unique(conn);
    Nv     = length(uverts);
    % determine distance between vertices
    D = nan(Nv);
    for ii=1:Nv
        for jj=ii+1:Nv
            D(jj,ii) = norm(XYZv(uverts(ii),:)-XYZv(uverts(jj),:));
        end
    end
    % find vertices that are closer than a threshhold value
    [jind,iind] = find(D <= tol);
    
    
    for ii = 1:length(jind)
        % replace one of the close-by vertices with their mean
        XYZv(uverts(iind(ii)),:) = mean(XYZv([uverts(iind(ii)),uverts(jind(ii))],:));
        % and replace reference to the other vertex in the connectivity list 
        conn(conn == uverts(jind(ii))) = uverts(iind(ii));
    end
    % to identify connected regions, determine the graph
    G = graph(categorical(conn(:,1)),categorical(conn(:,2)));
    % identify connected components of the graph
    [bins,binsize] = conncomp(G);
    Ng = unique(bins);
    
    for ii = 1:length(Ng)
        sG = subgraph(G,bins == ii);
        beg = find(degree(sG)==1,1,'first');
        
        if isempty(beg)
            beg = 1;
        end
        sorted_node_indices = dfsearch(sG,beg);
        sorted_indices = str2num(char(sG.Nodes{sorted_node_indices,1}));
            
        % convert to basis of plane   
        if (~nh(1)) & ~(nh(2)) & (nh(3)) % if the plane normal is parallel to the z-axis
            XYZp = [[0;0;sign(nh(3))],[1;0;0],[0;1;0]]^-1*(XYZv(sorted_indices,:)'-q);
        else 
            b1 = nh;
            b2 = cross([0;0;1],b1);
            b2 = b2/norm(b2);
            b3 = cross(b1,b2);
            XYZp = [b1,b2,b3]^-1*(XYZv(sorted_indices,:)'-q);
        end
        XYp{ii} = XYZp(2:3,:);
    end
end