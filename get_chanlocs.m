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

%% handle optional inputs
opts = cell2struct(varargin(2:2:end),varargin(1:2:end),2);
if ~isfield(opts,'anonymizeFace')
    opts.anonymizeFace = 0; end
if ~isfield(opts,'chanLabels')
    opts.chanLabels = {EEG.chanlocs(:).labels}'; end
if ~isfield(opts,'createMontageTemplate')
    opts.createMontageTemplate = 0; end
if ~isfield(opts,'fiducials')
    opts.fiducials = [];
    %opts.fiducials = [ 3.0022  -65.0460  109.4490;  76.5070  -39.5725   12.4267; -71.8769  -38.4040    9.6297 ]; end
if ~isfield(opts,'deleteTxtOutput')
    opts.deleteTxtOutput = 1; end
if ~isfield(opts, 'grayTextures')
    opts.grayTextures = 0; end
if ~isfield(opts,'moveElecInwards')
    opts.moveElecInwards = 0; end
if ~isfield(opts,'scannerAppName')
    opts.scannerAppName = ''; end

%% check 3D model input path for [.obj, .jpg,.mtl]
if nargin < 1
    help get_chanlocs;
    return;
end

%% selecting OBJ file
warndlg2( [ 'First select the OBJ file containing the scanned location.' 10 ...
            'This file is provided by your 3-D scanner' ])
[objFile, objPath] = uigetfile('*.obj', 'Select OBJ file');
if isequal(objFile, 0)
    disp('User abord.')
    return;
end

% additional options that depend on obj path
if ~isfield(opts,'templatePath') %#ok<ALIGN>
    opts.templatePath = [objPath, filesep, '..'];
    opts.templateSearch = 1;  end
if ~isfield(opts,'templateSaveName') %#ok<ALIGN>
    opts.templateSaveName = fullfile(opts.templatePath, 'montageTemplate.mat');
elseif isempty(regexp(opts.templateSaveName, '.mat','once'))
    fprintf('Appending file extension ".mat" to templateSaveName\n');
    opts.templateSaveName = [opts.templateSaveName,'.txt']; end
if ~isfield(opts,'saveName') %#ok<ALIGN>
    opts.saveName = fullfile(objPath, 'get_chanlocs.sfp');
elseif isempty(regexp(opts.saveName, '.sfp','once'))
    fprintf('Appending file extension ".sfp" to saveName\n');
    opts.saveName = [opts.saveName,'.sfp']; end

%% selecting template
choice = questdlg2(['It is easier to scan electrodes if you have already done so in another subject, and want to scan them again in no particular order.' 10 ...
    'The function will then assess electrode labels based on their position. Do you have a template you want to use?'], 'Template', ...
                           'No, I will select location in order', 'Yes, use MNI 10-20 reference locations', 'Yes, let me pick my .mat template file', 'No, I will select location in order');
if isempty(choice)
    disp('User abord.')
    return
end
opts.createMontageTemplate = 0;
switch choice
    case 'No, I will select location in order'
        opts.createMontageTemplate = 1;
    case 'Yes, use MNI 10-20 reference locations'
        opts.createMontageTemplate = 0;
        montageTemplate.refLocs = writeTemplateFromMNI(EEG);
    otherwise
        [fileNameTemplate, filePathTemplate] = uigetfile('*.mat', 'Select template file');
        if isequal(fileNameTemplate, 0)
            disp('User abord.')
            return
        end
        fprintf('Loading montage template from %s...\n', fullfile(filePathTemplate, fileNameTemplate));
        load('-mat', fullfile(filePathTemplate, fileNameTemplate))
end

%% anonymize face
if opts.anonymizeFace
    fprintf('Anonymizing face by replacing skintones with grey...\n')
    anonFace([objPath, filesep, getfield(dir([objPath filesep '*.jpg']),'name')]);
end

%% load model
fprintf('Loading 3D model in mm scale...\n')
head_surface = ft_read_headshape([objPath, filesep, getfield(dir([objPath filesep '*.obj']),'name')],...
    'unit', 'mm', 'grayTextures', opts.grayTextures);

%% locate fiducials and align
if isempty(opts.fiducials)
    fiducials = placeFiducials(head_surface);
else
    fiducials.elecpos = opts.fiducials;
end

fprintf('Using fiducials to align to CTF (mm) coordinates...\n')
cfg = []; cfg.method = 'fiducial'; cfg.coordsys = 'ctf';
cfg.fiducial.nas    = fiducials.elecpos(1,:); %position of NAS
cfg.fiducial.lpa    = fiducials.elecpos(2,:); %position of LHJ
cfg.fiducial.rpa    = fiducials.elecpos(3,:); %position of RHJ
[h, head_surface] = ft_meshrealign(cfg,head_surface);
rotatedFiducials = [fiducials.elecpos, ones(3,1)]*h';
rotatedFiducials(:,end) = [];

%% location selection
cfg = [];
cfg.channel = opts.chanLabels;

if opts.createMontageTemplate == 1
    warndlg2('Select electrode locations in the order electrodes are defined (use "r" to cancel selection)')
    elec = placeTemplateElectrodes(cfg,head_surface);
    elec.elecpos = fixupNewTemplate(cfg.channel, elec.elecpos, head_surface);
else
    cfg.montageTemplate = montageTemplate; cfg.refLocs = montageTemplate.refLocs; %#ok<NODEF>
    warndlg2('Select electrode locations in no particular order (use "r" to cancel selection)')
    elec = placeElectrodes(cfg,head_surface);
end

%% autoMapElectrodes
if ~opts.createMontageTemplate
    [elec.elecpos, affineTransformedRefLocs] = autoMapElectrodes(montageTemplate.refLocs, elec.elecpos);
    elec.elecpos = fixupMapping(montageTemplate.refLocs, elec.elecpos, affineTransformedRefLocs, head_surface);
end

%% save montageTemplate.mat
if opts.createMontageTemplate
    fprintf('Saving electrodes to %s so they can be reused for another subject\n', opts.templateSaveName)
    montageTemplate = head_surface; montageTemplate.refLocs = elec.elecpos(:,:);
    save(opts.templateSaveName, 'montageTemplate');
end

%% move electrodes in towards scalp
if opts.moveElecInwards ~= 0
    try
        disp('Moving electrodes toward the scalp')
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

fprintf('Saving electrode locations to txt file %s for manual editing/checking...\n', opts.saveName)
fileID = fopen(opts.saveName,'w');
for ii = 1:length(opts.chanLabels)
    fprintf(fileID, '%6s %9.4f %9.4f %9.4f\n', opts.chanLabels{ii}, -elec.elecpos(ii,2), elec.elecpos(ii,1), elec.elecpos(ii,3));
    % fprintf(fileID, '%6s %9.4f %9.4f %9.4f\n', opts.chanLabels{ii}, elec.elecpos(ii,1), elec.elecpos(ii,2), elec.elecpos(ii,3));
end
fclose(fileID);

fprintf('Inserting locations into the current dataset...\n')
EEG.chanlocs = readlocs(opts.saveName,'filetype', 'sfp'); %'format',{'labels','X','Y','Z'});
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
