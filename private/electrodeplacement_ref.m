function [elec] = ft_electrodeplacement(cfg, varargin)

% FT_ELECTRODEPLACEMENT allows manual placement of electrodes on a MRI scan, CT scan
% or on a triangulated surface of the head. This function supports different methods.
%
% VOLUME - Navigate an orthographic display of a volume (e.g. CT or MRI scan), and
% assign an electrode label to the current crosshair location by clicking on a label
% in the eletrode list. You can undo the selection by clicking on the same label
% again. The electrode labels shown in the list can be prespecified using cfg.channel
% when calling ft_electrodeplacement. The zoom slider allows zooming in at the
% location of the crosshair. The intensity sliders allow thresholding the image's low
% and high values. The magnet feature transports the crosshair to the nearest peak
% intensity voxel, within a certain voxel radius of the selected location. The labels
% feature displays the labels of the selected electrodes within the orthoplot. The
% global feature allows toggling the view between all and near-crosshair
% markers. The scan feature allows toggling between scans when another scan
% is given as input.
%
% HEADSHAPE - Navigate a triangulated scalp (for EEG) or brain (for ECoG) surface,
% and assign an electrode location by clicking on the surface. The electrode is
% placed on the triangulation itself.
%
% 1020 - Starting from a triangulated scalp surface and the nasion, inion, left and
% right pre-auricular points, this automatically constructs and follows contours over
% the surface according to the 5% system. Electrodes are placed at certain relative
% distances along these countours. This is an extension of the 10-20 standard
% electrode placement system and includes the 20%, 10% and 5% locations. See
% "Oostenveld R, Praamstra P. The five percent electrode system for high-resolution
% EEG and ERP measurements. Clin Neurophysiol. 2001 Apr;112(4):713-9" for details.
%
% Use as
%   [elec] = ft_electrodeplacement(cfg, ct)
%   [elec] = ft_electrodeplacement(cfg, ct, mri)
% where the input mri should be an anatomical CT or MRI volume, or
%   [elec] = ft_electrodeplacement(cfg, headshape)
% where the input headshape should be a surface triangulation.
%
% The configuration can contain the following options
%   cfg.method         = string representing the method for placing the electrodes
%                        'mri'             interactively locate electrodes in a MRI or CT scan
%                        'headshape'       interactively locate electrodes on a head surface
%                        '1020'            automatically place electrodes on a head surface
%
% The following options apply to the mri method
%   cfg.parameter      = string, field in data (default = 'anatomy' if present in data)
%   cfg.channel        = Nx1 cell-array with selection of channels (default = {'1' '2' ...})
%   cfg.elec           = struct containing previously placed electrodes (this overwrites cfg.channel)
%   cfg.clim           = color range of the data (default = [0 1], i.e. the full range)
%   cfg.magtype        = string representing the 'magnet' type used for placing the electrodes
%                        'peakweighted'    place electrodes at weighted peak intensity voxel (default)
%                        'troughweighted'  place electrodes at weighted trough intensity voxel
%                        'peak'            place electrodes at peak intensity voxel (default)
%                        'trough'          place electrodes at trough intensity voxel
%                        'weighted'        place electrodes at center-of-mass
%   cfg.magradius      = number representing the radius for the cfg.magtype based search (default = 3)
%
% The following options apply to the 1020 method
%   cfg.fiducial.nas   = 1x3 vector with coordinates
%   cfg.fiducial.ini   = 1x3 vector with coordinates
%   cfg.fiducial.lpa   = 1x3 vector with coordinates
%   cfg.fiducial.rpa   = 1x3 vector with coordinates
%   cfg.feedback       = string, can be 'yes' or 'no' for detailled feedback (default = 'yes')
%
% See also FT_ELECTRODEREALIGN, FT_VOLUMEREALIGN, FT_VOLUMESEGMENT, FT_PREPARE_MESH

% Copyright (C) 2015-2017, Arjen Stolk, Sandon Griffin & Robert Oostenveld
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

% give the user instructions
disp('Use the mouse to click on the indicated scalp positions');
disp('Press "r" to remove the last point added');
disp('Press "+/-" to zoom in/out');
disp('Press "w/a/s/d" to rotate');
disp('Press "q" to quit and advance when done');
% open a figure
figure('Name','getchanlocs')

% plot the faces of the 2D or 3D triangulation
headshape = fixpos(varargin{1});
ft_plot_mesh(headshape);

% plot the reference model on new axis and link the views
axes('Position', [0.6 0.5 0.5 0.5]); axis off
ft_plot_mesh(cfg.montageTemplate)

modelFig = findobj( 'Type', 'Figure', 'Name', 'getchanlocs' );
figAxes = get(modelFig,'Children');
Link = linkprop([figAxes(1), figAxes(2)], 'view');
setappdata(gcf, 'viewLink',Link);

% rotate3d on
xyz = ft_select_point3d(headshape, cfg.channel, cfg.refLocs, 'multiple', true, 'marker', '*');
numelec = size(xyz, 1);

% construct the output electrode structure
elec = keepfields(headshape, {'unit', 'coordsys'});
elec.elecpos = xyz;
for i=1:numelec
    elec.label{i} = cfg.channel{i,1};
end

function [selected] = ft_select_point3d(bnd, chanLabels, refLocs, varargin)

% FT_SELECT_POINT3D helper function for selecting one or multiple points on a 3D mesh
% using the mouse. It returns a list of the [x y z] coordinates of the selected
% points.
%
% Use as
%   [selected] = ft_select_point3d(bnd, ...)
%
% Optional input arguments should come in key-value pairs and can include
%   'multiple'    = true/false, make multiple selections, pressing "q" on the keyboard finalizes the selection (default = false)
%   'nearest'     = true/false (default = true)
%   'marker'      = character or empty, for example '.', 'o' or 'x' (default = [])
%   'markersize'  = scalar, the size of the marker (default = 10)
%   'markercolor' = character, for example 'r', 'b' or 'g' (default = 'k')
%
% Example
%   [pos, tri] = icosahedron162;
%   bnd.pos = pos;
%   bnd.tri = tri;
%   ft_plot_mesh(bnd)
%   camlight
%   ... do something here
%
% See also FT_SELECT_BOX, FT_SELECT_CHANNEL, FT_SELECT_POINT, FT_SELECT_RANGE, FT_SELECT_VOXEL

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

% get optional input arguments
marker      = ft_getopt(varargin, 'marker', []);
markersize  = ft_getopt(varargin, 'markersize', 10);
markercolor = ft_getopt(varargin, 'markercolor', 'k');

tmp = get(gcf,'Children'); refAxes = tmp(1); modelAxes = tmp(2); clear('tmp')
set(refAxes,'NextPlot','add'); set(modelAxes,'NextPlot','add');

selected = zeros(0,3);
done = false;
az = 0; el = 0;
view(az,el);

[sphereX, sphereY, sphereZ] = sphere;

ref_hs = surf(refAxes, sphereX*10 + refLocs(size(selected,1)+1,1),...
    sphereY*10 +refLocs(size(selected,1)+1,2),...
    sphereZ*10 +refLocs(size(selected,1)+1,3));
set(ref_hs, 'LineStyle', 'none', 'FaceColor',[1 0 0], 'FaceAlpha', 0.5);

supAxes = axes('pos',[0 0.95 1 1],'visible','off');
supText = text(supAxes,.5,0,['Select ' char(chanLabels(1))],...
    'FontSize',get(gcf,'defaultaxesfontsize')+4,...
    'horizontalalignment','center');


while ~done
    k = waitforbuttonpress;
    p = select3d(modelAxes);
    if k == 1 %checks if waitforbuttonpress was a key
        key = get(gcf,'CurrentCharacter'); % which key was pressed (if any)?
        if strcmp(key, 'q')
            % finished selecting points
            if size(selected,1)==size(chanLabels,1)
                done = true;
                fprintf('Location selection finished.\n')
            else
                warning('"q" press detected but points are missing!')
            end
        elseif strcmp(key, 'r')
            % remove last point
            if ~isempty(selected)
                if ~isempty(marker)
                    delete(findobj('marker', '*'));
                    hs = plot3(modelAxes, selected(1:end-1,1), selected(1:end-1,2), selected(1:end-1,3), [markercolor marker]);
                    set(hs, 'MarkerSize', markersize);
                end
                fprintf('Removed %s at [%9.4f %9.4f %9.4f]\n', char(chanLabels(size(selected,1))),...
                    selected(end,1), selected(end,2), selected(end,3));
                selected = selected(1:end-1,:);

                delete(ref_hs)
                ref_hs = surf(refAxes, sphereX*10 + refLocs(size(selected,1)+1,1),...
                    sphereY*10 +refLocs(size(selected,1)+1,2),...
                    sphereZ*10 +refLocs(size(selected,1)+1,3));
                set(ref_hs, 'LineStyle', 'none', 'FaceColor',[1 0 0], 'FaceAlpha', 0.5);
                set(supText,'String', ['Select ', char(chanLabels(size(selected,1)+1))]);
            end
        elseif strcmp(key,'+')
            zoom(modelAxes, 1.1)
        elseif strcmp(key,'-')
            zoom(modelAxes, 0.9)
        elseif strcmp(key,'w')
            el = el-6;
            view(modelAxes, az,el)
        elseif strcmp(key,'a')
            az = az+6;
            view(modelAxes, az,el)
        elseif strcmp(key,'s')
            el = el+6;
            view(modelAxes, az,el)
        elseif strcmp(key,'d')
            az = az-6;
            view(modelAxes, az,el)
        end
    else   % a new point was selected
        if size(selected,1)+1>size(chanLabels,1)
            fprintf(['Number of points selected exceed number of channels!\n'...
                'Last selected point was not added.\n'...
                'Press "q" to finish or "r" to remove last added point.\n'])
        elseif isempty(p)
            fprintf('Click not registered. Make sure points selected are within boundaries of 3D Model\n')
        elseif size(selected,1)+1==size(chanLabels,1)
            selected(end+1,:) = p;
            cla(refAxes);
            set(supText,'String', 'Press "q" to advance');
            fprintf('Selected %s at [%9.4f %9.4f %9.4f] (%d/%d) channels \n', char(chanLabels(size(selected,1))),...
                selected(end,1), selected(end,2), selected(end,3), size(selected,1),size(chanLabels,1));
            fprintf(['All channels now have locations.\n',...
                'Press "q" to quit and advance or "r" to remove last added point.\n'])
        else
            selected(end+1,:) = p;
            set(supText,'String', ['Select ', char(chanLabels(size(selected,1)+1))]);
            fprintf('Selected %s at [%9.4f %9.4f %9.4f] (%d/%d) channels \n', char(chanLabels(size(selected,1))),...
                selected(end,1), selected(end,2), selected(end,3), size(selected,1),size(chanLabels,1));
            if ~isempty(marker)&&~isempty(p)
                hs = plot3(modelAxes, selected(end,1), selected(end,2), selected(end,3), [markercolor marker]);
                set(hs, 'MarkerSize', markersize);
                
                delete(ref_hs)
                ref_hs = surf(refAxes, sphereX*10 + refLocs(size(selected,1)+1,1),...
                    sphereY*10 +refLocs(size(selected,1)+1,2),...
                    sphereZ*10 +refLocs(size(selected,1)+1,3));
                set(ref_hs, 'LineStyle', 'none', 'FaceColor',[1 0 0], 'FaceAlpha', 0.5);
            end
        end
    end
end

function [pout, vout, viout, facevout, faceiout]  = select3d(obj)
%SELECT3D(H) Determines the selected point in 3-D data space.
%  P = SELECT3D determines the point, P, in data space corresponding
%  to the current selection position. P is a point on the first
%  patch or surface face intersected along the selection ray. If no
%  face is encountered along the selection ray, P returns empty.
%
%  P = SELECT3D(H) constrains selection to graphics handle H and,
%  if applicable, any of its children. H can be a figure, axes,
%  patch, or surface object.
%
%  [P V] = SELECT3D(...), V is the closest face or line vertex
%  selected based on the figure's current object.
%
%  [P V VI] = SELECT3D(...), VI is the index into the object's
%  x,y,zdata properties corresponding to V, the closest face vertex
%  selected.
%
%  [P V VI FACEV] = SELECT3D(...), FACE is an array of vertices
%  corresponding to the face polygon containing P and V.
%
%  [P V VI FACEV FACEI] = SELECT3D(...), FACEI is the row index into
%  the object's face array corresponding to FACE. For patch
%  objects, the face array can be obtained by doing
%  get(mypatch,'faces'). For surface objects, the face array
%  can be obtained from the output of SURF2PATCH (see
%  SURF2PATCH for more information).
%
%  RESTRICTIONS:
%  SELECT3D supports surface, patch, or line object primitives. For surface
%  and patches, the algorithm assumes non-self-intersecting planar faces.
%  For line objects, the algorithm always returns P as empty, and V will
%  be the closest vertex relative to the selection point.
%
%  Example:
%
%  h = surf(peaks);
%  zoom(10);
%  disp('Click anywhere on the surface, then hit return')
%  pause
%  [p v vi face facei] = select3d;
%  marker1 = line('xdata',p(1),'ydata',p(2),'zdata',p(3),'marker','o',...
%                 'erasemode','xor','markerfacecolor','k');
%  marker2 = line('xdata',v(1),'ydata',v(2),'zdata',v(3),'marker','o',...
%                 'erasemode','xor','markerfacecolor','k');
%  marker2 = line('erasemode','xor','xdata',face(1,:),'ydata',face(2,:),...
%                 'zdata',face(3,:),'linewidth',10);
%  disp(sprintf('\nYou clicked at\nX: %.2f\nY: %.2f\nZ: %.2f',p(1),p(2),p(3)'))
%  disp(sprintf('\nThe nearest vertex is\nX: %.2f\nY: %.2f\nZ: %.2f',v(1),v(2),v(3)'))
%
%  Version 1.2 2-15-02
%  Copyright Joe Conti 2002
%  Send comments to jconti@mathworks.com
%
%  See also GINPUT, GCO.

% Output variables
pout = [];
vout = [];
viout = [];
facevout = [];
faceiout = [];

% other variables
ERRMSG = 'Input argument must be a valid graphics handle';
isline = false;
isperspective = false;

% Parse input arguments
if nargin<1
    obj = gco;
end

if isempty(obj) || ~ishandle(obj) || length(obj)~=1
    ft_error(ERRMSG);
end

% if obj is a figure
if strcmp(get(obj,'type'),'figure')
    fig = obj;
    ax = get(fig,'currentobject');
    currobj = get(fig,'currentobject');
    
    % bail out if not a child of the axes
    if ~strcmp(get(get(currobj,'parent'),'type'),'axes')
        return;
    end
    
    % if obj is an axes
elseif strcmp(get(obj,'type'),'axes')
    ax = obj;
    fig = get(ax,'parent');
    currobj = get(fig,'currentobject');
    currax = get(currobj,'parent');
    
    % Bail out if current object is under an unspecified axes
    if ~isequal(ax,currax)
        return;
    end
    
    % if obj is child of axes
elseif strcmp(get(get(obj,'parent'),'type'),'axes')
    currobj = obj;
    ax = get(obj,'parent');
    fig = get(ax,'parent');
    
    % Bail out
else
    return
end

axchild = currobj;
obj_type = get(axchild,'type');
is_perspective = strcmp(get(ax,'projection'),'perspective');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get projection transformation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% syntax not supported in old versions of MATLAB
[a b] = view(ax);
xform = viewmtx(a,b);
if is_perspective
    ft_warning('%s does not support perspective axes projection.',mfilename);
    d = norm(camtarget(ax)-campos(ax))
    P = [1 0 0 0;
        0 1 0 0;
        0 0 1 0;
        0 0 -1/d 1];
    xform = P*xform;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get vertex, face, and current point data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cp = get(ax,'currentpoint')';

% If surface object
if strcmp(obj_type,'surface')
    % Get surface face and vertices
    fv = surf2patch(axchild);
    vert = fv.vertices;
    faces = fv.faces;
    
    % If patch object
elseif strcmp(obj_type,'patch')
    vert = get(axchild,'vertices');
    faces = get(axchild,'faces');
    
    % If line object
elseif strcmp(obj_type,'line')
    xdata = get(axchild,'xdata');
    ydata = get(axchild,'ydata');
    zdata = get(axchild,'zdata');
    vert = [xdata', ydata',zdata'];
    faces = [];
    isline = true;
    
    % Ignore all other objects
else
    return;
end

% Add z if empty
if size(vert,2)==2
    vert(:,3) = zeros(size(vert(:,2)));
    if isline
        zdata = vert(:,3);
    end
end

% NaN and Inf check
nan_inf_test1 = isnan(faces) | isinf(faces);
nan_inf_test2 = isnan(vert) | isinf(vert);
if any(nan_inf_test1(:)) || any(nan_inf_test2(:))
    ft_warning('%s does not support NaNs or Infs in face/vertex data.',mfilename);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normalize for data aspect ratio %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dar = get(ax,'DataAspectRatio');

ncp(1,:) = cp(1,:)./dar(1);
ncp(2,:) = cp(2,:)./dar(2);
ncp(3,:) = cp(3,:)./dar(3);
ncp(4,:) = ones(size(ncp(3,:)));

nvert(:,1) = vert(:,1)./dar(1);
nvert(:,2) = vert(:,2)./dar(2);
nvert(:,3) = vert(:,3)./dar(3);
nvert(:,4) = ones(size(nvert(:,3)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transform data to view space %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xvert = xform*nvert';
xcp = xform*ncp;

if is_perspective % normalize 4th dimension
    xcp(1,:) = xcp(1,:)./xcp(4,:);
    xcp(2,:) = xcp(2,:)./xcp(4,:);
    xcp(3,:) = xcp(3,:)./xcp(4,:);
    xcp(4,:) = xcp(4,:)./xcp(4,:);
    
    xvert(1,:) = xvert(1,:)./xvert(4,:);
    xvert(2,:) = xvert(2,:)./xvert(4,:);
    xvert(3,:) = xvert(3,:)./xvert(4,:);
    xvert(4,:) = xvert(4,:)./xvert(4,:);
end

% Ignore 3rd & 4th dimensions for crossing test
xvert(4,:) = [];
xvert(3,:) = [];
xcp(4,:) = [];
xcp(3,:) = [];

% For debugging
% if 0
%     ax1 = getappdata(ax,'testselect3d');
%     if isempty(ax1) | ~ishandle(ax1)
%         fig = figure;
%         ax1 = axes;
%         axis(ax1,'equal');
%         setappdata(ax,'testselect3d',ax1);
%     end
%     cla(ax1);
%     patch('parent',ax1,'faces',faces,'vertices',xvert','facecolor','none','edgecolor','k');
%     line('parent',ax1,'xdata',xcp(1,2),'ydata',xcp(2,2),'zdata',0,'marker','o','markerfacecolor','r','erasemode','xor');
% end

% Translate vertices so that the selection point is at the origin.
xvert(1,:) = xvert(1,:) - xcp(1,2);
xvert(2,:) = xvert(2,:) - xcp(2,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simple algorithm (almost naive algorithm!) for line objects %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isline
    
    % Ignoring line width and marker attributes, find closest
    % vertex in 2-D view space.
    d = xvert(1,:).*xvert(1,:) + xvert(2,:).*xvert(2,:);
    [val i] = min(d);
    i = i(1); % enforce only one output
    
    % Assign output
    vout = [ xdata(i) ydata(i) zdata(i)];
    viout = i;
    
    return % Bail out early
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform 2-D crossing test (Jordan Curve Theorem) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find all vertices that have y components less than zero
vert_with_negative_y = zeros(size(faces));
face_y_vert = xvert(2,faces);
ind_vert_with_negative_y = find(face_y_vert<0);
vert_with_negative_y(ind_vert_with_negative_y) = true;

% Find all the line segments that span the x axis
is_line_segment_spanning_x = abs(diff([vert_with_negative_y, vert_with_negative_y(:,1)],1,2));

% Find all the faces that have line segments that span the x axis
ind_is_face_spanning_x = find(any(is_line_segment_spanning_x,2));

% Ignore data that doesn't span the x axis
candidate_faces = faces(ind_is_face_spanning_x,:);
vert_with_negative_y = vert_with_negative_y(ind_is_face_spanning_x,:);
is_line_segment_spanning_x = is_line_segment_spanning_x(ind_is_face_spanning_x,:);

% Create line segment arrays
pt1 = candidate_faces;
pt2 = [candidate_faces(:,2:end), candidate_faces(:,1)];

% Point 1
x1 = reshape(xvert(1,pt1),size(pt1));
y1 = reshape(xvert(2,pt1),size(pt1));

% Point 2
x2 = reshape(xvert(1,pt2),size(pt2));
y2 = reshape(xvert(2,pt2),size(pt2));

% Cross product of vector to origin with line segment
cross_product_test = -x1.*(y2-y1) > -y1.*(x2-x1);

% Find all line segments that cross the positive x axis
crossing_test = (cross_product_test==vert_with_negative_y) & is_line_segment_spanning_x;

% If the number of line segments is odd, then we intersected the polygon
s = sum(crossing_test,2);
s = mod(s,2);
ind_intersection_test = find(s~=0);

% Bail out early if no faces were hit
if isempty(ind_intersection_test)
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plane/ray intersection test %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform plane/ray intersection with the faces that passed
% the polygon intersection tests. Grab the only the first
% three vertices since that is all we need to define a plane).
% assuming planar polygons.
candidate_faces = candidate_faces(ind_intersection_test,1:3);
candidate_faces = reshape(candidate_faces',1,numel(candidate_faces));
vert = vert';
candidate_facev = vert(:,candidate_faces);
candidate_facev = reshape(candidate_facev,3,3,length(ind_intersection_test));

% Get three contiguous vertices along polygon
v1 = squeeze(candidate_facev(:,1,:));
v2 = squeeze(candidate_facev(:,2,:));
v3 = squeeze(candidate_facev(:,3,:));

% Get normal to face plane
vec1 = v2-v1;
vec2 = v3-v2;
crs = cross(vec1,vec2);
mag = sqrt(sum(crs.*crs));
nplane(1,:) = crs(1,:)./mag;
nplane(2,:) = crs(2,:)./mag;
nplane(3,:) = crs(3,:)./mag;

% Compute intersection between plane and ray
cp1 = cp(:,1);
cp2 = cp(:,2);
d = cp2-cp1;
dp = dot(-nplane,v1);

%A = dot(nplane,d);
A(1,:) = nplane(1,:).*d(1);
A(2,:) = nplane(2,:).*d(2);
A(3,:) = nplane(3,:).*d(3);
A = sum(A,1);

%B = dot(nplane,pt1)
B(1,:) = nplane(1,:).*cp1(1);
B(2,:) = nplane(2,:).*cp1(2);
B(3,:) = nplane(3,:).*cp1(3);
B = sum(B,1);

% Distance to intersection point
t = (-dp-B)./A;

% Find "best" distance (smallest)
[tbest ind_best] = min(t);

% Determine intersection point
pout = cp1 + tbest .* d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign additional output variables %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout>1
    
    % Get face index and vertices
    faceiout = ind_is_face_spanning_x(ind_intersection_test(ind_best));
    facevout = vert(:,faces(faceiout,:));
    
    % Determine index of closest face vertex intersected
    facexv = xvert(:,faces(faceiout,:));
    dist = sqrt(facexv(1,:).*facexv(1,:) +  facexv(2,:).*facexv(2,:));
    min_dist = min(dist);
    min_index = find(dist==min_dist);
    
    % Get closest vertex index and vertex
    viout = faces(faceiout,min_index);
    vout = vert(:,viout);
end


