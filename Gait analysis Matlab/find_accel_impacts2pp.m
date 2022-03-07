function [acc_pks,acc_pk_idx] = find_accel_impacts2pp(impacts,fsrTime,window,offset,accTime,accData,fsrData)
% returns the peak values and index of peaks for each impact
% acc_pks = [x peak, y peak, z peak, R = 1 or L = 0]
% acc_pk_idx = same as pks but index
% 11/22/21

Nimpacts = length(impacts(:,1));

acc_pks = zeros(Nimpacts,4); % cols: accX,accY,accZ peaks, whichfoot
acc_pk_idx = zeros(Nimpacts,4); % cols: accX,accY,accZ peak idx, whichfoot

for i=1:Nimpacts
    heel_start_t = fsrTime(impacts(i,1));
    window_start = heel_start_t - window*offset;
    window_end = window_start + window;
    starti = findTindex(window_start,accTime);
    endi = findTindex(window_end,accTime);
    if isempty(endi)
        endi = length(accTime);
    end
    footlabel = impacts(i,4);
    
    % accData has 12 cols
    accX = accData(starti:endi, (footlabel-1)*3 + 1);
    accY = accData(starti:endi, (footlabel-1)*3 + 2);
    accZ = accData(starti:endi, (footlabel-1)*3 + 3);
    [mx,ix] = max(abs(accX));
    [my,iy] = max(abs(accY));
    [mz,iz] = max(abs(accZ));
    acc_pks(i,:) = [mx,my,mz,footlabel];
    acc_pk_idx(i,:) = [ix+starti-1,iy+starti-1,iz+starti-1,footlabel];
    
end

% visualize
for i = 1:4 % 4 legs
    figure; subplot(2,1,1)
    plot(fsrTime, fsrData(:,i))
    hold on
    idx = find(impacts(:,4) == i);
    plot(fsrTime(impacts(idx,1)), fsrData(impacts(idx,1),i), 'rx')
    subplot(2,1,2)
    plot(accTime, abs(accData(:,(i-1)*3 + 2))) % y value only
    hold on
    idx = find(acc_pks(:,4)==i);
    plot(accTime(acc_pk_idx(idx,2)),acc_pks(idx,2),'rx')
    title("FSR and accel for foot " + i)
end
    
end