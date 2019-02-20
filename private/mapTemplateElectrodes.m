% map electrode positions from a working model to those from a template 
% model by minimizing euclidean distance (pdist2)
%
% inputs
% refLocs - reference electrode positions (from template model)
% scrambeldLocs - electrode positions (from working model)
% warningCost - minimized cost (total squred distance in mm) threshold
%               warning message printed if exceeded
% outputs:
% arrangedElecs - scrambledElec positions arranged in order of refLoc
% electrodes/ default channel ordering

function rearrangedLocs = mapTemplateElectrodes(refLocs, scrambledLocs, warningCost)
% munkres hungarian assignment
[initialReorderingIdx,T] = munkres(pdist2(refLocs,scrambledLocs));
rearrangedLocs = scrambledLocs(initialReorderingIdx,:);

flag_reorderingIdxUpdated = 1;
while flag_reorderingIdxUpdated == 1
    % affine transform
    % templateLocs = workingLocs*affineTransform + error %long Locs matrices
    templateLocs = [        refLocs ones(size(refLocs,1),1)];
    workingLocs  = [rearrangedLocs ones(size(refLocs,1),1)];
    A = workingLocs\templateLocs;
    
    affineTransformedLocs = workingLocs*A;
    affineTransformedLocs = affineTransformedLocs(:,1:3);
    [currentReorderingIdx,T] = munkres(pdist2(refLocs,affineTransformedLocs));
    if any(initialReorderingIdx ~= currentReorderingIdx)
        rearrangedLocs = rearrangedLocs(currentReorderingIdx,:);
        initialReorderingIdx = currentReorderingIdx;
    else
        flag_reorderingIdxUpdated = 0;
    end
end

if T >= warningCost
    warning(['Unpexectedly high distances found when automatically mapping '...
        'selected channels to template locations. Please double-check channel'...
        'labels manually or restart process and select each channel in order.'])
end
