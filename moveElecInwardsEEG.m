function EEG = moveElecInwardsEEG(EEG, electrodeHeight)
cfg.moveinward = electrodeHeight;
cfg.elec.coordsys = 'bti';
cfg.elec.unit = 'mm';
cfg.elec.elecpos = [[EEG.chanlocs.X]; [EEG.chanlocs.Y]; [EEG.chanlocs.Z]]';
elec = moveElecInwards(cfg);
EEG.chanlocs = keepfields(EEG.chanlocs, {'labels','type','urchan','ref'});

elecposCell = num2cell(elec.elecpos);
[EEG.chanlocs.X] = deal(elecposCell{:,1});
[EEG.chanlocs.Y] = deal(elecposCell{:,2});
[EEG.chanlocs.Z] = deal(elecposCell{:,3});