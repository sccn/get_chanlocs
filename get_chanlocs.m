% get_chanlocs(): Populate EEG.chanlocs using Wavefront .obj created by 3D-scanners.
% Currently tested to use models captured by the itSeez3D app and Structure scanner.
% For more information, visit https://sccn.ucsd.edu/wiki/get_chanlocs
% The github repository is located at https://github.com/cll008/get_chanlocs
% 
% FieldTrip toolbox functions are adapted to localize electrodes. See the original process at 
% http://www.fieldtriptoolbox.org/tutorial/electrode
%
% Usage:
%   >>EEG = get_chanlocs(EEG, objPath, 'key1', value1, ..., 'keyN', valueN);
%
% Inputs:
%   EEG     - EEGLAB EEG structure
%   objPath - full filepath to folder containing [*.obj, *.jpg, *.mtl] files
%
% Optional Inputs (value 0 or false to turn off):
%   'anonymizeFace'   - (Default = 0) Overwrite objPath .jpg file to replace 
%                        skintones with grey to anonymize subject's face
%   'chanLabels'      - (Default = {EEG.chanlocs(1,:).labels}') Label names
%                        for EEG channels (and misc sensors) to be localized.
%                        Default channel label list is extracted from EEG recording. 
%   'createMontageTemplate' - (Default = 0) Create and save rotated model with
%                        electrode label-location pairings to use as reference
%   'deleteTxtOutput' - (Default = 1) Delete text file containing electrode labels and 
%                        locations after importing to EEGLAB .set file.
%   'grayTextures'    - (Default = 0) Wrap gray textures around head model
%                        instead of textures prescribed in a separate .jpg
%                        file. Automatically set to 1 when no .jpg is found.
%   'moveElecInwards' - (Default = 0) Move electrode locations towards scalp
%                        (in mm) to adjust for cap and electrode well thickness.
%                        Negative numbers result in an outward move, away from the scalp.
%                        Scalar input moves all electrodes by the input value; if electrodes
%                        are offset by different heights, use a vector input with size 
%                        [1 x numberOfChannels]. e.g. [1 2.5 1 0] will move
%                        channels 1 and 3 by 1mm, channel 2 by 2.5mm, and leave channel 5 unmodified.
%   'templatePath'    - (Default = [objPath, filesep, '..']) Full filepath to montage template
%                       .mat file, or default to pop-up dialogue after searching model parent directory
%   'templateSaveName'- (Default = [templatePath, filesep, 'montageTemplate.mat'] filename (including path)
%                        of output file to save rotated model and locations to be used as template
%   'saveName'        - (Default = [objPath, filesep, 'get_chanlocs.txt']) Full filename 
%                       (including path) of output file to save electrode labels
%                        and locations. Imported into EEGLAB using readlocs(). Can be set to 
%                        automatically delete after import (see Optional Inputs: deleteTxtOutput)
%   'scannerAppName'  - (Default = []) Scanner app name (e.g. itSeez3D or Occipital Scanner)
%                        to be stored in EEG.chaninfo.get_chanlocs.scannerApp
%   'min_dim'         - (Default = 200) If the minimum range of the scan in
%                        any dimension is less than the MIN_DIM, the scan is unlikely to be 
%                        of the huamn head. Therefore, all dimensions will be scaled up by one order.
%                       To turn off set 'min_dim' 0.
%   'max_dim'         - (Default = 1000) If the maximum range of the scan in
%                        any dimension is greater than the MAX_DIM, the scan is unlikely to be 
%                        of the huamn head. Therefore, all dimensions will be scaled down by one order.
%                       To turn off set 'max_dim' 0.
%
% See also:
%   readlocs, ft_read_headshape, ft_electrodeplacement, ft_meshrealign
%   For FieldTrip (ft_*) functions, see http://www.fieldtriptoolbox.org
%   for the documentation and details.
%
% Author: 
%   Clement Lee, Swartz Center for Computational Neuroscience, 
%   Institute for Neural Computation, UC San Diego
%
% History: 
%   Github repo @ https://github.com/cll008/get_chanlocs
%   01 Feb 2018 v1.35 CL. Now supports reference model. Created repository
%   26 Jan 2018 v1.32 CL. try/catch on moveElecInwards for now so that process continues and locations are saved.
%   25 Jan 2018 v1.31 CL. Addressing issues with mex files when solid_angle.m (for moveElecInwards) fails. 
%   25 Jan 2018 v1.30 CL. Version check for string() when writing text file.
%   23 Jan 2018 v1.20 CL. Rename from GetChanLocs to getChanLocs. Default moveElecInwards 12.5-> 7.5.
%   19 Jan 2018	v1.10 CL. Switch fiducials to L/RHT (from L/RPA). Removed [] for channel / fiducial names. 
%   20 Dec 2017 v1.00 CL. Special thanks to M. Milham and L. Ai (@Child Mind Institute); N. Langer (@University of Zurich);
%					     and H. Tanaka (@Japan Advanced Institute of Science and Technology) for interest in testing.
%   4  Dec 2017 v0.10 Clement Lee. Created.
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

function EEG = get_chanlocs(EEG, objPath, varargin)
[~, versnum] = eeg_getversion;
if versnum < 2020
    error('Please use a newer version of EEGLAB (v2020.0+)')
end

%% check 3D model input path for [.obj, .jpg,.mtl]
if nargin < 1
	help get_chanlocs;
	return;
elseif mod(nargin,2) == 1
    help get_chanlocs;
    error('Please check input format.')
elseif ~(size(dir([objPath filesep '*.obj']),1) == 1&&...
         size(dir([objPath filesep '*.mtl']),1) == 1)
    error('Please input path to folder containing 1 set of [.obj and .mtl] files.')
end

%% handle optional inputs
opts = cell2struct(varargin(2:2:end),varargin(1:2:end),2);
if ~isfield(opts,'anonymizeFace')
    opts.anonymizeFace = 0; end
if ~isfield(opts,'chanLabels')
    opts.chanLabels = {EEG.chanlocs(:).labels}'; end
if ~isfield(opts,'createMontageTemplate')
    opts.createMontageTemplate = 0; end
if ~isfield(opts,'deleteTxtOutput')
    opts.deleteTxtOutput = 1; end
if ~isfield(opts, 'grayTextures')
    opts.grayTextures = 0; end
if ~isfield(opts, 'min_dim')
    opts.min_dim = 200; end
if ~isfield(opts, 'max_dim')
    opts.max_dim = 1000; end
if ~isfield(opts,'moveElecInwards')
    opts.moveElecInwards = 0; end
if ~isfield(opts,'templatePath') %#ok<ALIGN>
    opts.templatePath = [objPath, filesep, '..'];
    opts.templateSearch = 1;  end
if ~isfield(opts,'templateSaveName') %#ok<ALIGN>
    opts.templateSaveName = [opts.templatePath, filesep, 'montageTemplate.mat'];
elseif isempty(regexp(opts.templateSaveName, '.mat','once'))
    fprintf('Appending file extension ".mat" to templateSaveName\n');
    opts.templateSaveName = [opts.templateSaveName,'.txt']; end
if ~isfield(opts,'saveName') %#ok<ALIGN>
    opts.saveName = [objPath, filesep, 'get_chanlocs.txt'];
elseif isempty(regexp(opts.saveName, '.txt','once'))
    fprintf('Appending file extension ".txt" to saveName\n');
    opts.saveName = [opts.saveName,'.txt']; end
if ~isfield(opts,'scannerAppName')
    opts.scannerAppName = ''; end

%% anonymize face
if opts.anonymizeFace
fprintf('Anonymizing face by replacing skintones with grey...\n')
anonFace([objPath, filesep, getfield(dir([objPath filesep '*.jpg']),'name')]);
end

%% load model
fprintf('Loading 3D model in mm scale...\n')
head_surface = ft_read_headshape([objPath, filesep, getfield(dir([objPath filesep '*.obj']),'name')],...
    'unit', 'mm','grayTextures', opts.grayTextures);

% While ft_convert_units is meant to infer the units if they are not
% explicitly provided, the implementation is so unintuituve that I would
% rather import every scan as mm, and then scale them here to mm.

SCALE_UP = 0;
if min(range(head_surface.pos,1)) < opts.min_dim
    warning("the scan range is less than the minimum accepatble range for human head (200mm), scaling up the scan by one order")
    SCALE_UP = 1;
    head_surface.pos = head_surface.pos * 10;
end
if max(range(head_surface.pos,1)) > opts.max_dim
    if SCALE_UP
        error("The scan is heavily out of proportion, one dim is an order of magnitude smaller than the other, please redo the scan")
    end
    warning("the scan range is greater than the maximum accepatble range for human head (1000mm), scaling down the scan by one order")
    head_surface.pos = head_surface.pos / 10;
end

%% locate fiducials and align
fiducials = placeFiducials(head_surface);

fprintf('Using fiducials to align to CTF (mm) coordinates...\n')
cfg = []; cfg.method = 'fiducial'; cfg.coordsys = 'ctf';
cfg.fiducial.nas    = fiducials.elecpos(1,:); %position of NAS
cfg.fiducial.lpa    = fiducials.elecpos(2,:); %position of LHJ
cfg.fiducial.rpa    = fiducials.elecpos(3,:); %position of RHJ
[h, head_surface] = ft_meshrealign(cfg,head_surface);
rotatedFiducials = [fiducials.elecpos, ones(3,1)]*h';
rotatedFiducials(:,end) = [];

%% load reference template montage
if ~isfield(opts, 'templateSearch')
    load(opts.templatePath);
else
    refMats = dir([opts.templatePath, filesep, '*.mat']); 
    if isempty(refMats)
        choice = no_template;
        switch choice
            case 'Create new template'
                opts.createMontageTemplate = 1;
            case 'Use MNI template'
                opts.createMontageTemplate = 0;
                montageTemplate.refLocs = writeTemplateFromMNI(EEG);
%             case 'No'
%                 opts.createMontageTemplate = 0;
        end
    elseif size(refMats,1) == 1
        choice = one_template;
        switch choice
            case 'Load saved template'
                fprintf('Loading montage template from %s...\n', [opts.templatePath, filesep, refMats.name]);
                load([opts.templatePath, filesep, refMats.name])
            case 'Create new template'
                opts.createMontageTemplate = 1;
        end
    elseif size(refMats,1) > 1
        choice = multi_template(refMats);
        if strcmp(choice,'Create new template')
            opts.createMontageTemplate = 1;
        else
            load([opts.templatePath, filesep, choice]);
        end
    end
end
clear refMats choice

%% location selection
cfg = [];
cfg.channel = opts.chanLabels;

if opts.createMontageTemplate == 1
    fprintf('Select electrode locations for the new montage template...\n')
    elec = placeTemplateElectrodes(cfg,head_surface);
    elec.elecpos = fixupNewTemplate(cfg.channel, elec.elecpos, head_surface);
else
    cfg.montageTemplate = montageTemplate; cfg.refLocs = montageTemplate.refLocs; %#ok<NODEF>
    fprintf('Select electrode locations...\n')
    elec = placeElectrodes(cfg,head_surface);
end

%% autoMapElectrodes
if ~opts.createMontageTemplate
    [elec.elecpos, affineTransformedRefLocs] = autoMapElectrodes(montageTemplate.refLocs, elec.elecpos);
    elec.elecpos = fixupMapping(montageTemplate.refLocs, elec.elecpos, affineTransformedRefLocs, head_surface);
end

%% save montageTemplate.mat
if opts.createMontageTemplate
    fprintf('Saving reference .mat to %s\n', opts.templateSaveName)
    montageTemplate = head_surface; montageTemplate.refLocs = elec.elecpos(:,:);
    save(opts.templateSaveName, 'montageTemplate');
end

%% move electrodes in towards scalp
if opts.moveElecInwards ~= 0
    try
        cfg = [];
        cfg.elec = elec;
        cfg.moveinward = opts.moveElecInwards;
        elec = moveElecInwards(cfg);
    catch e
        fprintf(e.message)
        warning('Failed to move electrode positions. Saving unmodified coordinates...')
    end
end

%% format (labels','X','Y','Z') and save ascii for import. delete file afterwards if requested
% append rotated fiducials
opts.chanLabels(end+1:end+3,:) = {'nas','lhj','rhj'};
elec.elecpos(end+1:end+3,:) = rotatedFiducials;

fprintf('Writing electrode locations to txt file...\n')
fileID = fopen(opts.saveName,'w');
v = ver('Matlab');
if str2double(v.Version)>=9.1
    fprintf(fileID,'%6s %9.4f %9.4f %9.4f\n', [string(opts.chanLabels') ; elec.elecpos(:,:)']);
else
    for ii = 1:length(opts.chanLabels)
    fprintf(fileID, '%6s %9.4f %9.4f %9.4f\n', opts.chanLabels{ii}, elec.elecpos(ii,:));
    end
end
fclose(fileID);

fprintf('Importing locations with readlocs()...\n')
EEG.chanlocs = readlocs(opts.saveName,'format',{'labels','X','Y','Z'});
EEG = eeg_checkchanlocs(EEG);

if opts.deleteTxtOutput == 1
    fprintf('Deleting .txt file with electrode locations...\n')
    delete(opts.saveName)
end

% add get_chanlocs info to EEG.chaninfo
pathstr = fileparts(which('get_chanlocs'));
EEG.chaninfo.get_chanlocs.version = pathstr(end-3:end);
EEG.chaninfo.get_chanlocs.scannerApp = opts.scannerAppName;
EEG.chaninfo.get_chanlocs.transformMatrix = h;

fprintf('Electrode localization by get_chanlocs finished!\n')
end
