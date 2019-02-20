function elec = moveElecInwards(cfg)
% Move electrode locations towards scalp. Core function ft_electroderealign
% used to move electrodes towards scalp to account for cap and electrode
% well thickness. Augmented for the get_chanlocs EEGLAB plugin to receive
% non-scalar distance input (for caps with electrodes at different heights)
% Scalar input moves all electrodes by the input value; vector input with size 
% [1 x numberOfChannels]. e.g. [1 2.5 1 0] will move channels 1 and 3 by 1mm,
% channel 2 by 2.5mm, and leave channel 5 unmodified.
%
% Function handles created to accommodate direct calling from EEGLAB menu
% to shrink electrode locations without get_chanlocs process.
% 
% 12/11/2018 
% Clement Lee, Swartz Center for Computational Neuroscience,
% Institute for Neural Computation, UC San Diego
elec = cfg.elec;
cfg.method = 'moveinward';
cfg.elec.label = cell(1,size(elec.elecpos,1));

if size(cfg.moveinward,2) == 1
    cfg.moveinward = cfg.moveinward*ones(size(elec.elecpos,1),1)'; end
if size(cfg.moveinward,2) ~= size(elec.elecpos,1)
    error('moveElecInwards input vector length does not match number of electrodes'); end

[uniqueHeights, ~, subsetIdx] = unique(cfg.moveinward);
for ii = 1:length(uniqueHeights) % take care of one subset at a time
    cfg.moveinward = uniqueHeights(ii);
    tmpElec = ft_electroderealign(cfg);
    fprintf('Moving electrode set %d towards scalp by %2.2f mm...\n', ii, cfg.moveinward)
    elec.elecpos(subsetIdx == ii,:) = tmpElec.elecpos(subsetIdx == ii,:);
end
