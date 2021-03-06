function [elec] = placeElectrodes(cfg, head_surface)

% Based on ft_electrodeplacement(Copyright 2015-2017) by Arjen Stolk, Sandon Griffin & Robert Oostenveld
%
% Modified by:
%   Clement Lee, Swartz Center for Computational Neuroscience, 
%   Institute for Neural Computation, UC San Diego, 2018
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

% give the user instructions
disp('Use the mouse to select positions');
disp('Press "r" to remove the last point added');
disp('Press "+/-" to zoom in/out');
disp('Press "w/a/s/d" to rotate');
disp('Press "q" to quit and advance when done');
% open a figure
figure('Name','getchanlocs')

% plot
headshape = fixpos(head_surface);
ft_plot_mesh(headshape);
hold on
if isfield(cfg,'fixupLocs')
    [sphereX, sphereY, sphereZ] = sphere;
    for locIdx = 1:size(cfg.fixupLocs,1)
        fixup_hs = surf(sphereX*7.5 + cfg.fixupLocs(locIdx,1),...
            sphereY*7.5 + cfg.fixupLocs(locIdx,2),...
            sphereZ*7.5 + cfg.fixupLocs(locIdx,3));
        set(fixup_hs, 'LineStyle', 'none', 'FaceColor',[1 0 0], 'FaceAlpha', 0.5);
    end

else
    cfg.fixupLocs = zeros(0,3); 
end

% select
xyz = ft_select_point3d(cfg);
close gcf

% construct the output electrode structure
elec = keepfields(headshape, {'unit', 'coordsys'});
elec.elecpos = xyz;


function [selected] = ft_select_point3d(cfg)

% FT_SELECT_POINT3D helper function for selecting one or multiple points on a 3D mesh
% using the mouse. It returns a list of the [x y z] coordinates of the selected
% points.
% See also FT_SELECT_BOX, FT_SELECT_CHANNEL, FT_SELECT_POINT, FT_SELECT_RANGE, FT_SELECT_VOXEL

% get optional input arguments
marker = '*'; markersize = 10; markercolor = 'k';
refLocs = cfg.refLocs;
selected = cfg.fixupLocs;

modelAxes = get(gcf,'Children');
supAxes = axes('pos',[0 0.95 1 1],'visible','off');
supText = text(supAxes,.5,0, ['Point selection progress: ',...
    num2str(size(selected,1)) ' out of ' num2str(size(refLocs,1))],...
    'FontSize',get(gcf,'defaultaxesfontsize')+4, 'horizontalalignment','center');
set(modelAxes,'NextPlot','add'); 

done = false;
az = 0; el = 0;
view(az,el);

while ~done
    k = waitforbuttonpress;
    p = select3d(modelAxes);
    if k == 1 %checks if waitforbuttonpress was a key
        key = get(gcf,'CurrentCharacter'); % which key was pressed (if any)?
        if strcmp(key, 'q')
            % finished selecting points
            if size(selected,1)==size(refLocs,1)
                done = true;
                fprintf('Location selection finished.\n')
            else
                warning(['"q" press detected but number of points selected'...
                    'is not equal to number of channels in the template!'])
            end
        elseif strcmp(key, 'r')
            % remove last point
            if ~isempty(selected)
                if ~isempty(marker)
                    delete(findobj('marker', '*'));
                    hs = plot3(modelAxes, selected(1:end-1,1), selected(1:end-1,2),...
                        selected(1:end-1,3), [markercolor marker]);
                    set(hs, 'MarkerSize', markersize);
                end
                fprintf('Removed point at [%9.4f %9.4f %9.4f]\n',...
                    selected(end,1), selected(end,2), selected(end,3));
                selected = selected(1:end-1,:);
                set(supText,'String', ['Point selection progress: ',...
                    num2str(size(selected,1)) ' out of ' num2str(size(refLocs,1))]);
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
        if size(selected,1)+1 > size(refLocs,1)
            fprintf(['Number of points selected exceed number of channels in the template!\n'...
                'Last selected point was not added.\n'...
                'Press "q" to finish or "r" to remove last added point.\n'])
        elseif isempty(p)
            fprintf('Click not registered. Make sure points selected are within boundaries of 3D Model\n')
        else
            selected(end+1,:) = p;
            fprintf('Selected point at [%9.4f %9.4f %9.4f] (%d/%d) channels \n', ...
                selected(end,1), selected(end,2), selected(end,3), size(selected,1), size(refLocs,1));
            
            if ~isempty(marker)&&~isempty(p)
                hs = plot3(modelAxes, selected(end,1), selected(end,2), selected(end,3), [markercolor marker]);
                set(hs, 'MarkerSize', markersize);
            end
            
            if size(selected,1) == size(refLocs,1)
                set(supText,'String', 'Press "q" to quit and advance');
                fprintf(['Number of points selected now matches the number of channels in the template.\n',...
                    'Press "q" to quit and advance or "r" to remove last added point.\n'])
            else
                set(supText,'String', ['Point selection progress: ',...
                    num2str(size(selected,1)) ' out of ' num2str(size(refLocs,1))]);
            end
        end
    end
end