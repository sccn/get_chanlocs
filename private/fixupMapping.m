function rearrangedLocs = fixupMapping(refLocs, rearrangedLocs, affineTransformedRefLocs, head_surface)

figure; plotElectrodePairings(affineTransformedRefLocs, rearrangedLocs)

%interactive utility
disp('Press "r" to remove a selected location');
disp('Press "s" to select new location(s)');
disp('Press "c" to compute new assignment and show updated plot');
disp('Press "q" to quit and advance');
done = false;

while ~done
    k = waitforbuttonpress;
    if k == 1
        key = get(gcf,'CurrentCharacter');
        if strcmp(key, 'q')
            if size(rearrangedLocs,1) == size(refLocs,1)
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
            if size(refLocs,1) - size(rearrangedLocs,1) < 1 
                fprintf(['Number of channels in current model will exceed that in the template model\n'...
                    'Please remove a channel before adding any more']);
            else
                cfg.refLocs = refLocs;
                cfg.fixupLocs = rearrangedLocs;
                elec = placeElectrodes(cfg, head_surface);
                rearrangedLocs = elec.elecpos(:,:);
                disp('Press "c" to compute new assignment and show updated plot');
            end
        elseif strcmp(key,'c')
            if size(rearrangedLocs,1) == size(refLocs,1)
                [rearrangedLocs, affineTransformedRefLocs] = autoMapElectrodes(refLocs, rearrangedLocs);
                clf; plotElectrodePairings(affineTransformedRefLocs, rearrangedLocs);
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

function plotElectrodePairings(refLocs, rearrangedLocs)
%% visual confirmation plot to check validity of auto mapping
shape1 = refLocs;
shape2 = rearrangedLocs;

plot3(shape1(:,1),shape1(:,2),shape1(:,3),'o')
hold on; plot3(shape2(:,1), shape2(:,2), shape2(:,3),'ko')
for n = 1:size(shape1,1)
    plot3([shape1(n,1), shape2(n,1)],...
        [shape1(n,2), shape2(n,2)],...
        [shape1(n,3), shape2(n,3)],'r','linewidth',1.5)
    text(shape1(n,1)+2.5,shape1(n,2)-2.5,shape1(n,3)+2.5, num2str(n),'HorizontalAlignment','center',...
        'VerticalAlignment','middle','Color',[0 0 0])
end
axis equal vis3d off
suptitle('Channel Assignment Result'); view(135,35)
legend('Transformed Template Locations','Current Model Locations','Assigned Pairing and Channel Index','Location','northeast')
