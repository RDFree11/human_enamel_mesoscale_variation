function [QryPts] = diffpol2cart_getQueryPts(twotheta,chi,offset,tilt,Sdd)
%pol2cart - Converts a polar 2D diffraction pattern image to its cartesian
%form in 2-theta space, taking inputs:
%   'Interpolant' - Pass in an interpolant of the Diffraction pattern to be
%   transformed (created using griddedInterpolant) on raw image.
%   'twotheta' - 1-D array containing twotheta values for transform
%   (implicitly defines resolution of interpolation along twotheta)
%   'chi' - 1-D array containing chi values for transform (implicitly defines resolution of interpolation along chi)
%   'offset' - [dxc,dyx] array containing the offset of the beamcenter from
%   that used to construct the interpolant in the detector coordinate
%   system
%   indexed from the bottom left
%   'tilt' - [tilt rotation, tilt angle] array containing the rotation of
%   the tilt plane in degrees from the ideal +x axis and the angle of
%   rotation out of the ideal detector plane along the tilt plane.
%   'Sdd' - sample-to-detector distance in pixels.


% define grid of rays in terms of Bragg angle and azimuthal angle
[Chiq,TwoThetaq] = ndgrid(chi,twotheta);

% convert query rays into unit vectors in Cartesian coordinates (beam center of detector is
% (0,0,0) and beam direction is [0,0,-1]
Rq(1,:,:) = sind(TwoThetaq) .* cosd(Chiq); 
Rq(2,:,:) = sind(TwoThetaq) .* sind(Chiq); 
Rq(3,:,:) = cosd(TwoThetaq);

% Compute orientation of detector plane
k = [sind(tilt(1));-cosd(tilt(1));0]; %vector definition of rotation axis to transform from ideal detector to actual detector normal.
K=[0 0 k(2); 0 0 -k(1); -k(2) k(1) 0]; %K matrix for Rodriguez formula, rotating about the k axis
I=[1 0 0; 0 1 0; 0 0 1]; %identity
R = I + sind(tilt(2))*K + (1-cosd(tilt(2)))*K^2; %Rotation matrix from Rodriguez formula
n_det = R*[0;0;1]; %normal vector of real detector plane

% Determine coordinates of intersection of query rays with detector plane
N = repmat(n_det,1,numel(chi),numel(twotheta)); %for efficient calculation of dot products, create array of normal vectors
NU = (Sdd*n_det(3))./squeeze(dot(N,Rq,1)); %compute scaling factor. (0,0,-Sdd) + nu*Rq_hat = (x,y,z) of detector intersection in lab coordinates
Lq(1,:,:) = NU.*squeeze(Rq(1,:,:));
Lq(2,:,:) = NU.*squeeze(Rq(2,:,:));
Lq(3,:,:) = NU.*squeeze(Rq(3,:,:))-Sdd;

% Transform coordinates into the basis of the detector (relative to the
% beam center)
QryPts = (R^-1)*reshape(Lq,3,[]);

% Adjust query points based on offset of the beam center from the one assumed by the
% interpolant (needs to be done after basis transformation b/c beamshift
% corresponds to movement of origin of coordinate system for previous
% vector math
QryPts(1,:) = QryPts(1,:,:) + offset(1); %x offset
QryPts(2,:) = QryPts(2,:,:) + offset(2); %y offset
end