function [rearrangedLocs, affineTransformedRefLocs, currentReorderingIdx] = autoMapElectrodes(refLocs, scrambledLocs)
% map selected electrodes from a working model to those from a template 
% model by minimizing euclidean distance (pdist2). Uses affine mapping to
% transform template locations and munkres assignment. scrambledLocs
% may be incomplete (less than total number of electrodes), in which case 
% rearrangedLocs will fill in missing electrodes with affine transformed
% template locations.
%
% inputs
% refLocs - reference electrode positions (from template model)
% scrambeldLocs - electrode positions (from working model) 
%
% outputs:
% arrangedElecs - scrambledElec positions arranged in order of refLoc
%                 electrodes (default channel ordering)

% munkres-hungarian assignment
[initialReorderingIdx,~] = munkres(pdist2(scrambledLocs, refLocs));
tmpRefLocs = refLocs(initialReorderingIdx,:);

flag_reorderingIdxUpdated = 1;
while flag_reorderingIdxUpdated == 1
    % affine transform
    % solve affine transform based on available information
    workingLocs  = [scrambledLocs ones(size(scrambledLocs,1),1)];    
    templateLocs = [   tmpRefLocs ones(size(scrambledLocs,1),1)];
    A = templateLocs\workingLocs;
    
    % apply transform to full & correctly-ordered reference. re-munkres
    tmpRefLocs = [refLocs ones(size(refLocs,1),1)]*A;
    tmpRefLocs = tmpRefLocs(:,1:3);
    [currentReorderingIdx,~] = munkres(pdist2(scrambledLocs, tmpRefLocs));
    
    if any(initialReorderingIdx ~= currentReorderingIdx)
        tmpRefLocs = refLocs(currentReorderingIdx,:);
        initialReorderingIdx = currentReorderingIdx;
    else
        flag_reorderingIdxUpdated = 0;
        affineTransformedRefLocs = tmpRefLocs;
        rearrangedLocs = affineTransformedRefLocs;
        rearrangedLocs(currentReorderingIdx,:) = scrambledLocs;
    end
end

%% test results
% using 111 ratings from 3 users on 13 head models (SCCN 2018 rTMS basic reaserch project)
% note for developer: see incompleteAffineTest.m for details

% When assigning selected locations from each rating to MNI template:
%      0     0     0     0     0     0     0     0     0     0
%      0     0     0     0     0     0     0     0     1     2
%      2     3     1     0     3     3     6     5    10    10
%      12    12    16    20    28    28    32    42    51    59
%     103   103   109   162   177   244   299   340   404   463 
%     544   594   683   714   810   724   690   666   506   336   185
%
% 0 out of 5550 failures with 61 of 61 electrodes...
% 1  (0.02%) @ 42 / 61 electrodes...
% 10 (0.18%) @ 32 / 61 electrodes...
% 103 (1.9%) @ 11 / 61 electrodes...
% max 810 (14.73%) @ 6 / 61 electrodes...
% 
% When comparing choose 2 from 111 models (6105 comparisons total):
% 0           2           6          12          22          42          76         115         157       196
% 275         322         367         417         478         489         530         584      644         669
% 724         741         777         812         844         895         946         949         982        1021
% 1041        1097        1164        1173        1201        1290        1288        1374        1396        1457
% 1565        1560        1641        1718        1767        1885        1974        2061        2191        2278
% 2449        2541        2610        2671        2723        2707        2574        2475        2125        1718        1074
% 
% 0.02% error at 57 / 61 electrodes, 0.2% at 53 / 61 electrodes, 2% at 17 / 61 electrodes
