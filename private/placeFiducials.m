function [elec] = placeFiducials(head_surface)

% Based on ft_electrodeplacement (Copyright 2015-2017) by Arjen Stolk, Sandon Griffin & Robert Oostenveld
%
% Modified by:
%   Clement Lee, Swartz Center for Computational Neuroscience, 
%   Institute for Neural Computation, UC San Diego, 2019
%
% This original file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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
cfg.channel = {'NAS','LHJ','RHJ'}';

% give the user instructions
fprintf('Select Nasion, Left Helix/Tragus Junction, and Right Helix/Tragus Junction...\n')
disp('Use the mouse to click on the indicated fiducials');
disp('Press "r" to remove the last point added');
disp('Press "+/-" to zoom in/out');
disp('Press "w/a/s/d" to rotate');
disp('Press "q" to quit and advance when done');
% open a figure
figure('Name','getchanlocs')

% plot and select
headshape = fixpos(head_surface);
ft_plot_mesh(headshape);
xyz = ft_select_point3d(cfg);
close gcf

xyz_old = xyz;
figure('Name','getchanlocs')
ft_plot_mesh(headshape);
fprintf('Select fiducials again\n')
xyz = ft_select_point3d(cfg);
close gcf
pass = 0;

while ~pass
    if any(diag(pdist2(xyz_old,xyz)) > 2)
        xyz_old = xyz;
        figure('Name','getchanlocs')
        ft_plot_mesh(headshape);
        warning('Fiducial distance tolerance (2mm) exceeded. Select fiducials again.')
        xyz = ft_select_point3d(cfg);
        close gcf
    else
        fprintf('Fiducial distances within tolerance. Averaging...\n')
        xyz = (xyz+xyz_old)/2;
        pass = 1;
    end
end

% construct the output electrode structure
elec = keepfields(headshape, {'unit', 'coordsys'});
elec.elecpos = xyz;

function [selected] = ft_select_point3d(cfg)

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
marker = '*'; markersize = 10; markercolor = 'k';
chanLabels = cfg.channel;
modelAxes = get(gcf, 'Children');
set(modelAxes,'NextPlot','add');

done = false;
az = 0;
el = 0;
view(az,el);

selected = zeros(0,3);
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
                    hs = plot3(selected(1:end-1,1), selected(1:end-1,2), selected(1:end-1,3), [markercolor marker]);
                    set(hs, 'MarkerSize', markersize);
                end
                fprintf('Removed %s at [%9.4f %9.4f %9.4f]\n', char(chanLabels(size(selected,1))),...
                    selected(end,1), selected(end,2), selected(end,3));
                selected = selected(1:end-1,:);
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
    else        % a new point was selected
        if size(selected,1)+1>size(chanLabels,1)
            fprintf(['Number of points selected exceed number of fiducials!\n'...
                'Last selected point was not added.\n'...
                'Press "q" to finish or "r" to remove last added point.\n'])
        elseif isempty(p)
            fprintf('Click not registered. Make sure points selected are within boundaries of 3D Model\n')
        elseif size(selected,1)+1==size(chanLabels,1)
            selected(end+1,:) = p;
            set(supText,'String', 'Press "q" to advance');
            fprintf('Selected %s at [%9.4f %9.4f %9.4f] \n', char(chanLabels(size(selected,1))),...
                selected(end,1), selected(end,2), selected(end,3));
            fprintf(['All fiducials now have locations.\n',...
                'Press "q" to quit and advance or "r" to remove last added point.\n'])
        else
            selected(end+1,:) = p;
            set(supText,'String', ['Select ', char(chanLabels(size(selected,1)+1))]);
            fprintf('Selected %s at [%9.4f %9.4f %9.4f] \n', char(chanLabels(size(selected,1))),...
                selected(end,1), selected(end,2), selected(end,3));
        end
        if ~isempty(marker)&&~isempty(p)
            hs = plot3(selected(end,1), selected(end,2), selected(end,3), [markercolor marker]);
            set(hs, 'MarkerSize', markersize);
        end
    end
end