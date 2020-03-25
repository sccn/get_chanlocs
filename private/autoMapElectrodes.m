function [rearrangedLocs, affineTransformedLocs] = autoMapElectrodes(refLocs, scrambledLocs)
% map electrode positions from a working model to those from a template 
% model by minimizing euclidean distance (pdist2)
%
% inputs
% refLocs - reference electrode positions (from template model)
% scrambeldLocs - electrode positions (from working model)
%
% outputs:
% arrangedElecs - scrambledElec positions arranged in order of refLoc
%                 electrodes (default channel ordering)

% munkres-hungarian assignment
[initialReorderingIdx,~] = munkres(pdist2(refLocs,scrambledLocs));
rearrangedLocs = scrambledLocs(initialReorderingIdx,:);

flag_reorderingIdxUpdated = 1;
while flag_reorderingIdxUpdated == 1
    % affine transform
    % templateLocs = workingLocs*affineTransform + error %long Locs matrices
    templateLocs = [       refLocs ones(size(refLocs,1),1)];
    workingLocs  = [rearrangedLocs ones(size(refLocs,1),1)];
    A = workingLocs\templateLocs;
    
    affineTransformedLocs = workingLocs*A;
    affineTransformedLocs = affineTransformedLocs(:,1:3);
    [currentReorderingIdx,~] = munkres(pdist2(refLocs,affineTransformedLocs));
    if any(initialReorderingIdx ~= currentReorderingIdx)
        rearrangedLocs = rearrangedLocs(currentReorderingIdx,:);
        initialReorderingIdx = currentReorderingIdx;
    else
        flag_reorderingIdxUpdated = 0;
    end
end