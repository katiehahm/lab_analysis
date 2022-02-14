function [impacts] = findimpacts_fsr_accel2pp(fsrTime, mocapT, allmocap, fsrData,dist,min_thresh)
% finds the heel impacts from fsr data that includes accel
% by looking at right and left heel strikes, findpeaks, and filtering out
% wide peaks and peaks with no =0 between them
% impacts = [heel start idx, heel pk idx, heel pk mag, 1(R) or 0(L)]
% 9/15/21

impacts = [0, 0, 0, 0]; % initialize array

for i = 1:4
    [pks, locs, ~,~] = findpeaks(fsrData(:,i), 'MinPeakProminence',8,'Annotate','extents','MinPeakDistance',dist(i));
    
    % filter out peaks that don't have heel go close to 0 in between
    omit = [];
    for j = 2:length(pks)
        btw_pk = fsrData(locs(j-1):locs(j),i);
        if ~any(btw_pk < min_thresh(i))
            omit(end+1) = j;
        end
    end
    locs(omit) = [];
    pks(omit) = [];
    
    % find actual start of heel strike instead of peak
    N = length(locs);
    for j = 2:N
        window = fsrData(locs(j-1):locs(j),i);
        [~,idx] = min(window);
        window = movmean(window, 20);
        start_idx = aic_pick(window(idx:end),'to_peak') + idx + locs(j-1);
        impacts(end+1,:) = [start_idx, locs(j), pks(j), i];
    end
end

impacts(1,:) = []; % delete first row from initialization
impacts = sortrows(impacts,1);

% eliminate any impacts off the floor/turning points
x_min = -3500;
x_max = 3600;
elim_i = [1]; % first step is always on the resonance

for i = 1:length(impacts)
    impact_time = fsrTime(impacts(i,1));
    mocap_i = findTindex(impact_time,mocapT);
    person_i = impacts(i,4);
    coord = allmocap(mocap_i, 3, person_i); % 3 is for z dir
    if coord > x_max || coord < x_min
        elim_i(end+1) = i;
    end
end

impacts(elim_i,:) = [];

% % visually check
% % fsr
% RorL = length(impacts(1,:)); % always the last column
% figure;
% subplot(2,1,1)
% hold on
% plot(fsrData(:,Mfsr('Lheel')))
% lefts = find(impacts(:,RorL) == 0);
% plot(impacts(lefts,1),fsrData(impacts(lefts,1),Mfsr('Lheel')),'rx','MarkerSize',12)
% plot(impacts(lefts,2),fsrData(impacts(lefts,2),Mfsr('Lheel')),'bx','MarkerSize',12)
% title('Left heel')
% 
% subplot(2,1,2)
% hold on
% plot(fsrData(:,Mfsr('Rheel')))
% rights = find(impacts(:,RorL) == 1);
% plot(impacts(rights,1),fsrData(impacts(rights,1),Mfsr('Rheel')),'rx','MarkerSize',8)
% plot(impacts(rights,2),fsrData(impacts(rights,2),Mfsr('Rheel')),'bx','MarkerSize',12)
% title('Right heel')
% 
% % pcb
% figure;
% plot(pcbTime,filt_pcbD(:,1))
% hold on
% plot(fsrTime(impacts(:,1)),0,'r.','MarkerSize',12)

% visually check
for i = 1:4
    figure;
    idxs = find(impacts(:,4) == i);
    hold on
    plot(fsrData(:,i))
    plot(impacts(idxs,1),fsrData(impacts(idxs,1),i), 'rx', 'MarkerSize', 12)
    plot(impacts(idxs,2),impacts(idxs,3),'bx','MarkerSize',12)
    title(['Foot ', string(i)])
end