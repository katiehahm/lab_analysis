function [acc_pks,acc_pk_idx] = find_accel_impacts(impacts,fsrTime,window,offset,accTime,accData,fsrData,Mfsr)
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
    RorLfoot = impacts(i,4);
    if RorLfoot == 1 % right 
        accX = accData(starti:endi,(Mfsr('Rheel')-1)*3 + 1);
        accY = accData(starti:endi,(Mfsr('Rheel')-1)*3 + 2);
        accZ = accData(starti:endi,(Mfsr('Rheel')-1)*3 + 3);
        [mx,ix] = max(abs(accX));
        [my,iy] = max(abs(accY));
        [mz,iz] = max(abs(accZ));
        acc_pks(i,:) = [mx,my,mz,RorLfoot];
        acc_pk_idx(i,:) = [ix+starti-1,iy+starti-1,iz+starti-1,RorLfoot];
    else
        accX = accData(starti:endi,(Mfsr('Lheel')-1)*3 + 1);
        accY = accData(starti:endi,(Mfsr('Lheel')-1)*3 + 2);
        accZ = accData(starti:endi,(Mfsr('Lheel')-1)*3 + 3);
        [mx,ix] = max(abs(accX));
        [my,iy] = max(abs(accY));
        [mz,iz] = max(abs(accZ));
        acc_pks(i,:) = [mx,my,mz,RorLfoot];
        acc_pk_idx(i,:) = [ix+starti-1,iy+starti-1,iz+starti-1,RorLfoot];
    end
end

% visualize
figure; subplot(4,1,1)
plot(fsrTime,fsrData(:,1))
hold on
idx = find(impacts(:,4)==0);
plot(fsrTime(impacts(idx,1)),0,'rx')
subplot(4,1,2)
plot(accTime, abs(accData(:,1)))
hold on
idx = find(acc_pks(:,4)==0);
plot(accTime(acc_pk_idx(idx,1)),acc_pks(idx,1),'rx')
subplot(4,1,3)
plot(accTime, abs(accData(:,2)))
hold on
idx = find(acc_pks(:,4)==0);
plot(accTime(acc_pk_idx(idx,2)),acc_pks(idx,2),'rx')
subplot(4,1,4)
plot(accTime, abs(accData(:,3)))
hold on
idx = find(acc_pks(:,4)==0);
plot(accTime(acc_pk_idx(idx,3)),acc_pks(idx,3),'rx')
title('Left foot fsr and accel peaks')

figure; subplot(4,1,1)
plot(fsrTime,fsrData(:,2))
hold on
idx = find(impacts(:,4)==1);
plot(fsrTime(impacts(idx,1)),0,'rx')
subplot(4,1,2)
plot(accTime, abs(accData(:,4)))
hold on
idx = find(acc_pks(:,4)==1);
plot(accTime(acc_pk_idx(idx,1)),acc_pks(idx,1),'rx')
subplot(4,1,3)
plot(accTime, abs(accData(:,5)))
hold on
idx = find(acc_pks(:,4)==1);
plot(accTime(acc_pk_idx(idx,2)),acc_pks(idx,2),'rx')
subplot(4,1,4)
plot(accTime, abs(accData(:,6)))
hold on
idx = find(acc_pks(:,4)==1);
plot(accTime(acc_pk_idx(idx,3)),acc_pks(idx,3),'rx')
title('Right foot fsr and accel peaks')
end