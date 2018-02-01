function [elec_realigned] = ft_electroderealign(cfg, elec_original)

% FT_ELECTRODEREALIGN rotates, translates, scales and warps electrode positions. The
% default is to only rotate and translate, i.e. to do a rigid body transformation in
% which only the coordinate system is changed. With the right settings if can apply
% additional deformations to the input sensors (e.g. scale them to better fit the
% skin surface). The different methods are described in detail below.
%
% INTERACTIVE - You can display the skin surface together with the electrode or
% gradiometer positions, and manually (using the graphical user interface) adjust the
% rotation, translation and scaling parameters, so that the electrodes correspond
% with the skin.
%
% FIDUCIAL - You can apply a rigid body realignment based on three fiducial
% locations. After realigning, the fiducials in the input electrode set (typically
% nose, left and right ear) are along the same axes as the fiducials in the template
% electrode set.
%
% TEMPLATE - You can apply a spatial transformation/deformation that automatically
% minimizes the distance between the electrodes or gradiometers and a template or
% sensor array. The warping methods use a non-linear search to minimize the distance
% between the input sensor positions and the corresponding template sensors.
%
% HEADSHAPE - You can apply a spatial transformation/deformation that automatically
% minimizes the distance between the electrodes and the head surface. The warping
% methods use a non-linear search to minimize the distance between the input sensor
% positions and the projection of the electrodes on the head surface.
%
% PROJECT - This projects all electrodes to the nearest point on the
% head surface mesh.
%
% MOVEINWARD - This moves all electrodes inward according to their normals
%
% Use as
%   [elec_realigned] = ft_sensorrealign(cfg)
% with the electrode or gradiometer details in the configuration, or as
%   [elec_realigned] = ft_sensorrealign(cfg, elec_orig)
% with the electrode or gradiometer definition as 2nd input argument.
%
% The configuration can contain the following options
%   cfg.method         = string representing the method for aligning or placing the electrodes
%                        'interactive'     realign manually using a graphical user interface
%                        'fiducial'        realign using three fiducials (e.g. NAS, LPA and RPA)
%                        'template'        realign the electrodes to match a template set
%                        'headshape'       realign the electrodes to fit the head surface
%                        'project'         projects electrodes onto the head surface
%                        'moveinward'      moves electrodes inward along their normals
%   cfg.warp          = string describing the spatial transformation for the template and headshape methods
%                        'rigidbody'       apply a rigid-body warp (default)
%                        'globalrescale'   apply a rigid-body warp with global rescaling
%                        'traditional'     apply a rigid-body warp with individual axes rescaling
%                        'nonlin1'         apply a 1st order non-linear warp
%                        'nonlin2'         apply a 2nd order non-linear warp
%                        'nonlin3'         apply a 3rd order non-linear warp
%                        'nonlin4'         apply a 4th order non-linear warp
%                        'nonlin5'         apply a 5th order non-linear warp
%                        'dykstra2012'     non-linear wrap only for headshape method, useful for projecting ECoG onto cortex hull
%                        'fsaverage'       surface-based realignment with the freesurfer fsaverage brain
%   cfg.channel        = Nx1 cell-array with selection of channels (default = 'all'),
%                        see  FT_CHANNELSELECTION for details
%   cfg.keepchannel    = string, 'yes' or 'no' (default = 'no')
%   cfg.fiducial       = cell-array with the name of three fiducials used for
%                        realigning (default = {'nasion', 'lpa', 'rpa'})
%   cfg.casesensitive  = 'yes' or 'no', determines whether string comparisons
%                        between electrode labels are case sensitive (default = 'yes')
%   cfg.feedback       = 'yes' or 'no' (default = 'no')
%
% The electrode positions can be present in the 2nd input argument or can be specified as
%   cfg.elec          = structure with electrode positions, see FT_DATATYPE_SENS
%   cfg.elecfile      = name of file containing the electrode positions, see FT_READ_SENS
%
% If you want to realign the EEG electrodes using anatomical fiducials, you should
% specify the target location of the three fiducials, e.g.
%   cfg.target.pos(1,:) = [110 0 0]     % location of the nose
%   cfg.target.pos(2,:) = [0  90 0]     % location of the left ear
%   cfg.target.pos(3,:) = [0 -90 0]     % location of the right ear
%   cfg.target.label    = {'NAS', 'LPA', 'RPA'}
%
% If you want to align EEG electrodes to a single or multiple template electrode sets
% (which will be averaged), you should specify the template electrode sets either as
% electrode structures (i.e. when they are already read in memory) or their file
% names using
%   cfg.target          = single electrode set that serves as standard
% or
%   cfg.target{1..N}    = list of electrode sets that will be averaged
%
% If you want to align EEG electrodes to the head surface, you should specify the head surface as
%   cfg.headshape      = a filename containing headshape, a structure containing a
%                        single triangulated boundary, or a Nx3 matrix with surface
%                        points
%
% If you want to align ECoG electrodes to the pial surface, you first need to compute
% the cortex hull with FT_PREPARE_MESH. dykstra2012 uses algorithm described in
% Dykstra et al. (2012, Neuroimage) in which electrodes are projected onto pial
% surface while minimizing the displacement of the electrodes from original location
% and maintaining the grid shape. It relies on the optimization toolbox.
%   cfg.method         = 'headshape'
%   cfg.warp           = 'dykstra2012'
%   cfg.headshape      = a filename containing headshape, a structure containing a
%                        single triangulated boundary, or a Nx3 matrix with surface
%                        points
%   cfg.feedback       = 'yes' or 'no' (default), feedback of the iteration procedure
%
% If you want to move the electrodes inward, you should specify
%   cfg.moveinward     = number, the distance that the electrode should be moved
%                        inward (negative numbers result in an outward move)
%
% If you want to align ECoG electrodes to the freesurfer average brain, you should
% specify the path to your headshape (e.g., lh.pial), and ensure you have the
% corresponding registration file (e.g., lh.sphere.reg) in the same directory.
% Moreover, the path to the local freesurfer home is required. Note that, because the
% electrodes are being aligned to the fsaverage brain, the corresponding brain should
% be also used when plotting the data, i.e. use freesurfer/subjects/fsaverage/surf/lh.pial
% rather than surface_pial_left.mat.
%   cfg.method         = 'headshape'
%   cfg.warp           = 'fsaverage'
%   cfg.headshape      = string, filename containing headshape (e.g. <path to freesurfer/surf/lh.pial>)
%   cfg.fshome         = string, path to freesurfer
%
% See also FT_READ_SENS, FT_VOLUMEREALIGN, FT_INTERACTIVEREALIGN,
% FT_DETERMINE_COORDSYS, FT_PREPARE_MESH

% Copyright (C) 2005-2015, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% the interactive method uses a global variable to get the data from the figure when it is closed
global norm

% set the defaults
cfg.warp          = ft_getopt(cfg, 'warp', 'rigidbody');
cfg.channel       = ft_getopt(cfg, 'channel',  'all');
cfg.feedback      = ft_getopt(cfg, 'feedback', 'no');
cfg.casesensitive = ft_getopt(cfg, 'casesensitive', 'no');
cfg.headshape     = ft_getopt(cfg, 'headshape', []);     % for triangulated head surface, without labels
cfg.coordsys      = ft_getopt(cfg, 'coordsys');          % this allows for automatic template fiducial placement
    
% the data can be passed as input arguments or can be read from disk
hasdata = exist('elec_original', 'var');

% get the electrode definition that should be warped
% if ~hasdata
%   elec_original = ft_fetch_sens(cfg);
% else
%   % the input electrodes were specified as second input argument
%   % or read from cfg.inputfile
% end
elec_original = cfg.elec;

% % ensure that channel and electrode positions are the same
% assert(isequaln(elec_original.elecpos, elec_original.chanpos), 'this function requires same electrode and channel positions');

% remember the original electrode locations and labels and do all the work with a
% temporary copy, this involves channel selection and changing to lower case
elec = elec_original;
label_original = elec_original.label;

% start with an empty structure, this will be returned at the end
norm = [];

  norm.label = elec.label;
  norm.elecpos = moveinward(elec.elecpos, cfg.moveinward);

% apply the spatial transformation to all electrodes, and replace the
% electrode labels by their case-sensitive original values
    % nothing to be done
    elec_realigned = norm;
    elec_realigned.label = label_original;


% the coordinate system is in general not defined after transformation
if isfield(elec_realigned, 'coordsys')
  elec_realigned = rmfield(elec_realigned, 'coordsys');
end

% the coordinate system remains the same
elec_realigned.coordsys = elec_original.coordsys;

% channel positions are identical to the electrode positions (this was checked at the start)
elec_realigned.chanpos = elec_realigned.elecpos;

% copy over unit, chantype, and chanunit information in case this was not already done
elec_realigned.unit = elec_original.unit;
end

function [pos] = moveinward(pos, move)
%This functions moves 'pos' inward according to their normals by 'move'
%units
propos = elproj(pos); % projection to 2D
tri = delaunay(propos); %creates delaunay triangulation of 2D plane, which will be used for the the 3D case
nor = normals(pos,tri); %compute normals of surface
ori = surfaceorientation(pos, tri, nor);
if ori==1
  % the normals are outward oriented
elseif ori==-1
  % the normals are inward oriented
  nor = -nor;
else
  error('cannot determine the orientation of the vertex normals');
end
pos = pos-move*nor; % moves pos inwards according to their normals
end

function [nrm] = normals(pnt, tri, opt)

% NORMALS compute the surface normals of a triangular mesh
% for each triangle or for each vertex
%
% [nrm] = normals(pnt, tri, opt)
% where opt is either 'vertex' or 'triangle'

% Copyright (C) 2002-2007, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

npnt = size(pnt,1);
ntri = size(tri,1);

% shift to center
pnt(:,1) = pnt(:,1)-mean(pnt(:,1),1);
pnt(:,2) = pnt(:,2)-mean(pnt(:,2),1);
pnt(:,3) = pnt(:,3)-mean(pnt(:,3),1);

v2 = pnt(tri(:,2),:) - pnt(tri(:,1),:);
v3 = pnt(tri(:,3),:) - pnt(tri(:,1),:);
nrm_tri = cross(v2, v3);

% compute vertex normals
nrm_pnt = zeros(npnt, 3);
for i=1:ntri
nrm_pnt(tri(i,1),:) = nrm_pnt(tri(i,1),:) + nrm_tri(i,:);
nrm_pnt(tri(i,2),:) = nrm_pnt(tri(i,2),:) + nrm_tri(i,:);
nrm_pnt(tri(i,3),:) = nrm_pnt(tri(i,3),:) + nrm_tri(i,:);
end
% normalise the direction vectors to have length one
nrm = nrm_pnt ./ (sqrt(sum(nrm_pnt.^2, 2)) * ones(1,3));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fast cross product to replace the MATLAB standard version
function [c] = cross(a,b)
c = [a(:,2).*b(:,3)-a(:,3).*b(:,2) a(:,3).*b(:,1)-a(:,1).*b(:,3) a(:,1).*b(:,2)-a(:,2).*b(:,1)];
end

function val = surfaceorientation(pnt, tri, ori)

% SURFACEORIENTATION returns 1 if the triangulated surface is outward
% oriented, -1 if it is inward oriented and 0 if the orientation cannot be
% determined.
%
% Use as
%   surfaceorientation(pnt, tri)
% or
%   surfaceorientation(pnt, tri, ori)

% Copyright (C) 2007, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

pnt(:,1) = pnt(:,1)-mean(pnt(:,1),1);
pnt(:,2) = pnt(:,2)-mean(pnt(:,2),1);
pnt(:,3) = pnt(:,3)-mean(pnt(:,3),1);

% FIXME there is a bug in solid_angle resulting in negative values where they should be positive and vice versa 
if all(sign(sum(pnt .* ori, 2))==1)
  % the normals are outward oriented
  val = 1;
elseif all(sign(sum(pnt .* ori, 2))==-1)
  % the normals are inward oriented
  val = -1;
elseif abs(sum(solid_angle(pnt, tri))+4*pi)<1000*eps
  % the normals are outward oriented
  val = 1;
elseif abs(sum(solid_angle(pnt, tri))-4*pi)<1000*eps
  % the normals are inward oriented
  val = -1;
else
  % cannot determine
  val = 0;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that applies the homogeneous transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new] = apply(transform, old)
old(:,4) = 1;
new = old * transform';
new = new(:,1:3);
end
