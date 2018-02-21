function [refHeadModel, refLocs] = prepareRefModel(refPath)
fprintf('Loading reference model...\n')
refHeadModel = ft_read_headshape(strcat(refPath, filesep, 'Model.obj'), 'unit','mm');

fprintf('Reading reference locations...\n')
refLocFile = [refPath, filesep, 'getChanLocs.txt'];
array = loadtxt(refLocFile,'verbose','off','blankcell','off');
array(:,1) = []; refLocs = cell2mat(array);

fprintf('Aligning reference to BTi coordinates...\n')
cfg = []; cfg.method = 'fiducial'; cfg.coordsys = 'bti';
cfg.fiducial.nas    = refLocs(end-2,:); %position of NAS
cfg.fiducial.lpa    = refLocs(end-1,:); %position of LHT
cfg.fiducial.rpa    = refLocs(end  ,:); %position of RHT
refHeadModel = ft_meshrealign(cfg,refHeadModel);
end