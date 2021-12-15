function rearrangedLocs = fixupNewTemplate(chanLabels, inputLocs, head_surface)

rearrangedLocs = inputLocs;
figure; plotNewTemplate(chanLabels, inputLocs, head_surface)

%interactive utility
disp('Press "r" to remove a selected location');
disp('Press "s" to select new location(s)');
disp('Press "c" to compute and show updated plot');
disp('Press "q" to quit and advance');
done = false;

while ~done
    k = waitforbuttonpress;
    if k == 1
        key = get(gcf,'CurrentCharacter');
        if strcmp(key, 'q')
            if size(rearrangedLocs,1) == size(chanLabels,1)
                close
                done = true;
            else
                disp('Number of channels selected does not equal that in the template model')
                disp('Channel assignment is incomplete!')
            end
        elseif strcmp(key, 'r')
            rmIdx = input(['Which location(s) should be removed?\n'...
                'Please enter the channel index/indices (e.g. 1 or [1 4 7]):']);
            rearrangedLocs(rmIdx,:) = [];
            disp('Channel(s) removed! Press "s" to select new location(s)');
        elseif strcmp(key,'s')
            if size(chanLabels,1) - size(rearrangedLocs,1) < 1 
                fprintf(['Number of channels in current model will exceed that in the template model\n'...
                    'Please remove a channel before adding any more']);
            else
                cfg.refLocs = chanLabels;
                cfg.fixupLocs = rearrangedLocs;
                elec = placeElectrodes(cfg, head_surface);
                rearrangedLocs = elec.elecpos(:,:);
                disp('Press "c" compute and show updated plot');
            end
        elseif strcmp(key,'c')
            if size(rearrangedLocs,1) == size(chanLabels,1)
                clf; plotNewTemplate(chanLabels, rearrangedLocs, head_surface)
                disp('Press "r" to remove a selected location');
                disp('Press "s" to select new location(s)');
                disp('Press "q" to quit and advance');
            else
                disp('Number of channels selected does not equal that in the template model')
                disp('Press "r" to remove a selected location');
                disp('Press "s" to select new location(s)');
            end
        end
    end
end

function plotNewTemplate(chanLabels, inputLocs, head_surface)
%% visual confirmation plot to check validity of auto mapping
headshape = fixpos(head_surface);
ft_plot_mesh(headshape);
hold on

[sphereX, sphereY, sphereZ] = sphere;
for locIdx = 1:size(inputLocs,1)
    fixup_hs = surf(sphereX*7.5 + inputLocs(locIdx,1),...
        sphereY*7.5 + inputLocs(locIdx,2),...
        sphereZ*7.5 + inputLocs(locIdx,3));
    set(fixup_hs, 'LineStyle', 'none', 'FaceColor',[1 0 0], 'FaceAlpha', 0.5);
    
    text(inputLocs(locIdx,1),inputLocs(locIdx,2),inputLocs(locIdx,3)+10,...
        [num2str(locIdx) '. ' chanLabels{locIdx}], 'HorizontalAlignment','Center', 'Color',[0 0 0])
end

axis equal vis3d off; view(135,35)
supAxes = axes('pos',[0 0.95 1 1],'visible','off');
text(supAxes,.5,0,['New template result'],...
    'FontSize',get(gcf,'defaultaxesfontsize')+4,...
    'horizontalalignment','center');
