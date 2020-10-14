% eegplugin_getchanlocs()
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Notes:
%   To create a new plugin, simply create a file beginning with "eegplugin_"
%   and place it in your eeglab folder. It will then be automatically
%   detected by eeglab. See also this source code internal comments.
%   For eeglab to return errors and add the function's results to
%   the eeglab history, menu callback must be nested into "try" and
%   a "catch" strings. 
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA

function vers = eegplugin_getchanlocs(fig, trystrs, catchstrs)

vers = 'get_chanlocs 3.00';
if nargin < 3
    error('eegplugin_getchanlocs requires 3 arguments');
end;

% create menu
toolsMenu = findobj(fig, 'tag', 'tools');
uimenu(toolsMenu, 'label', 'Locate electrodes in a 3D head image',...
    'callback', 'gui_getchanlocs','separator','on');