% check_plugin_vers() - checks if the local plug-in has the same version 
% number as the remote listing hosted on the EEGLAB extension manager.
%
% Usage:
% >> [flag_diff_vers, remote_plugin_vers] = check_plugin_vers(pluginName, versNumChar, verbose);
%
% Inputs:
%       pluginName  - (string) name of plug-in to check
%       versNumChar - (int) number of characters used to denote version at 
%                     the end of the parent directory's name. (e.g. 
%                     versNumChar = 4 for ../pluginName2.38/pluginName)
%       verbose     - (boolean, default 1) prints output to Command Window
% 
% Outputs:
%       flag_diff_vers     - (boolean) true if version numbers are different
%       remote_plugin_vers - (string)  version number of remote listing 
%
% Author: 
%   Clement Lee, Swartz Center for Computational Neuroscience, 
%   Institute for Neural Computation, UC San Diego
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

function [flag_diff_vers, remote_plugin_vers] = check_plugin_vers(pluginName,versNumChar,verbose)

if nargin < 1
	help check_plugin_vers;
	return;
elseif  ~exist('verbose','var')
    verbose = 1;
end
    
localPluginPath = which(pluginName);
localPluginDir  = fileparts(localPluginPath);
localPluginVers = localPluginDir(end-(versNumChar-1):end);

remotePluginStruct = plugin_getweb('process');

remotePluginNames = {remotePluginStruct.name};
remotePluginIdx = find(strcmp(remotePluginNames,pluginName),1);

if isempty(remotePluginIdx)
    fprintf('Could not find the plug-in listing on the extension maanger\n')
    return
else
    remote_plugin_vers = remotePluginStruct(remotePluginIdx).version;
    flag_diff_vers = strcmp(remote_plugin_vers, localPluginVers);
end

if verbose == 1
    if flag_diff_vers == 0
        warning(['Version number of ' pluginName ' does not match the listing in the extension manager!'])
        fprintf('Check File->Manage EEGLAB extensions\n')
    else
        fprintf(['Version number of ' pluginName ' matches the listing in the extension manager!\n'])
    end
end

