%% 4/18/22 experiment 4 all processing, do it section by section
clear all
close all
filepath = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\';
% ################# change ########################
intervention = 'both_regular2'; % regular1 limp1 limp2 weight1 weight2 regular2
% #################################################

% load 3 datasets
load([filepath, intervention])
load([filepath, 'fsr_', intervention])
T = readtable([filepath, 'mocap_', intervention]);

% constants
sensorN = 6;
Fs_pcb = 12800;
MCkeySet = {'Lx1','Ly1','Lz1','Rx1','Ry1','Rz1','Lx2','Ly2','Lz2','Rx2','Ry2','Rz2'};
MCvalueSet = [3 4 5 6 7 8 9 10 11 12 13 14]; % order of Mocap left/right based on keyset
Mmocap = containers.Map(MCkeySet, MCvalueSet);

Fs_mocap = 240;
[mocapT, mocapL1, mocapR1, mocapL2, mocapR2] = convertMocap2pp(T, Mmocap);
allmocap = cat(3, mocapL1,mocapR1,mocapL2,mocapR2);

[accData, accTime, fsrData, fsrTime] = clip_accelfsr_fromMocap2pp_ver2(fsrData, accData, trigger_signal, Fs_fsr, Fs_acc, Fs_trigger);
[pcbData, pcbTime] = clip_pcb_fromMocap(datas, times, sensorN);

% check this is correct. These numbers should be nearly the same:
mocapT(end)
fsrTime(end)
pcbTime(end)
accTime(end)

%% delete start and end of data when walking on/off
figure; plot(mocapT, mocapL2(:,1),'b.')
figure; plot(pcbTime, pcbData)

%% cont.

start_time = 11.2; % change ##
startidx_pcb = findTindex(start_time,pcbTime);
startidx_fsr = findTindex(start_time,fsrTime);
startidx_acc = findTindex(start_time,accTime);
startidx_mocap = findTindex(start_time,mocapT);

fsrData = fsrData(startidx_fsr:end,:);
fsrTime = fsrTime - fsrTime(startidx_fsr);
fsrTime = fsrTime(startidx_fsr:end);
pcbTime = pcbTime - pcbTime(startidx_pcb);
pcbTime = pcbTime(startidx_pcb:end,:);
pcbData = pcbData(startidx_pcb:end,:);
mocapT = mocapT - mocapT(startidx_mocap);
mocapT = mocapT(startidx_mocap:end);
allmocap = allmocap(startidx_mocap:end,:,:);
accTime = accTime - accTime(startidx_acc);
accTime = accTime(startidx_acc:end);
accData = accData(startidx_acc:end,:);


last_time = 95.8 - start_time; % change ##
lastidx_pcb = findTindex(last_time,pcbTime);
lastidx_fsr = findTindex(last_time,fsrTime);
lastidx_acc = findTindex(last_time,accTime);
lastidx_mocap = findTindex(last_time,mocapT);

fsrData = fsrData(1:lastidx_fsr,:);
fsrTime = fsrTime(1:lastidx_fsr);
pcbTime = pcbTime(1:lastidx_pcb);
pcbData = pcbData(1:lastidx_pcb,:);
mocapT = mocapT(1:lastidx_mocap);
allmocap = allmocap(1:lastidx_mocap,:,:);
accTime = accTime(1:lastidx_acc);
accData = accData(1:lastidx_acc,:);

mocapT(end)
fsrTime(end)
pcbTime(end)
accTime(end)

figure; plot(mocapT, allmocap(:,1,3),'b.')
figure; plot(pcbTime, pcbData)

%% clean up foot drops error in sensor jenny second half 5/21/22 (limp2)
% first run next section. Then run "clean up foot drops", then run this
% section, then run next section again.

% old_fsrData = fsrData;
% % set any drops under 5 = prev value
% for i = 1:length(fsrData(:,4))
%     curr_val = fsrData(i,4);
%     if curr_val < 6.5
%         fsrData(i,4) = fsrData(i-1,4);
%     end
% end

%%
% clean and filter pcb data
% filt_pcbD = lpf_data(pcbData);

% finding footfalls based on fsr heel data
L_dist = 180; % min distance between heel peaks
R_dist = 180;
min_threshL = 25; % value between peaks threshold, if not lower than this then omit peak
min_threshR = 25;

% min distance between heel peaks
dist = [180, 180, 180, 180]; % L1, R1, L2, R2
min_thresh = [25, 25, 25, 25]; % L1, R1, L2, R2

% nan_idx = find(fsrData(:,4) < 5);
% fsrData(nan_idx,4) = NaN;

impacts = findimpacts_fsr_accel2pp(fsrTime, mocapT, allmocap, fsrData,dist,min_thresh);

%% to clean up foot drops in sensor 5/21/22 

% % set all value after a peak to be flat zero
% wrong_impacts = find(impacts(:,4) == 4); % foot 4
% fsrL = length(fsrData(:,4));
% for i = 1:length(wrong_impacts)
%     curr_peak_idx = impacts(wrong_impacts(i),2);
%     fsrData(min(curr_peak_idx + 5,fsrL):min(curr_peak_idx + 140,fsrL),4) = 8.5;
% end
% figure; plot(fsrData(:,4))
% then run above section again



%% fix small errors in impacts 11/5/21
heel_start_wrong = [1754,11094,13574,14544,15758]; % these need to be same length
heel_start_right = [1827,11131,13704,14524,15795];

heel_pk_wrong = []; % index, these need to be same length
heel_pk_right = [];

delete_pks = [4208]; % index of peaks to delete

impacts = manual_fix_fsr2pp(impacts,fsrData,heel_start_wrong,heel_start_right,heel_pk_wrong,heel_pk_right,delete_pks);
impacts = sortrows(impacts,1);
%% mocap extract
[coordinates] = findimpacts_mocap2pp(impacts,fsrTime,mocapT,allmocap,true);

%% extract accel data
% gets the peak of abs value of accelerometer values for x,y,z directions
% 4th col in variable 'acc_pks' indicates whether it is LorR leg.
window = 0.3; % seconds
offset = 0.25;

[acc_pks,acc_pk_idx] = find_accel_impacts2pp(impacts,fsrTime,window,offset,accTime,accData,fsrData);

%% saving data to matlab
processedfilepath = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\';
filename = [processedfilepath,intervention];

save(filename,'filename','pcbTime','datas','pcbData', ...
    'fsrTime','fsrData','impacts','Fs_pcb','Fs_fsr','Fs_acc','Fs_trigger', ...
    'acc_pks','acc_pk_idx','accTime','accData',...
    'mocapT','allmocap','coordinates','Fs_mocap','sensorN')
disp(append("Saved as ", filename))
close all
%% FOOTFALL DETECTION #############################################################################################

%% find walk edges and segments

% eliminate impacts that are off the floor
elim_idx = [];
for i = 1:length(impacts)
    curr_coord = coordinates(i,1); % x value
    if curr_coord > 3.5 | curr_coord < -3.5
        elim_idx(end+1) = i;
    end
end
impacts(elim_idx,:) = [];
coordinates(elim_idx,:) = [];
acc_pk_idx(elim_idx,:) = [];
acc_pks(elim_idx,:) = [];


% label walk edges and segments
impactN = length(impacts);
walk_edges = zeros(1,impactN);
walk_edges(1) = -1;
walk_edges(end) = 1;
ID_labels = zeros(1,impactN);

% ID labels
o_idx = find(impacts(:,4) == 1 | impacts(:,4) == 2);
x_idx = find(impacts(:,4) == 3 | impacts(:,4) == 4);
ID_labels(o_idx) = 1;
ID_labels(x_idx) = 2;

% find walk_edges
for i = 2:impactN
    if (impacts(i,1) - impacts(i-1,1))/Fs_fsr > 1.5 % start of new walking segment
        walk_edges(i-1) = 1;
        walk_edges(i) = -1;
    end
end

% find segments
segments = find(walk_edges == -1);
segments(end+1) = impactN+1;

% find average step time for person o
step_o1_diff = [];
step_o2_diff = [];
real_o_idx = find(impacts(:,4) == 1 | impacts(:,4) == 2);
real_o_times = impacts(real_o_idx,1)./Fs_fsr;
for x = 2:length(real_o_times)
    diff = real_o_times(x)-real_o_times(x-1);
    if diff < 1.5 % not a turning point
        if impacts(real_o_idx(x),4) == 1 % left foot
            step_o1_diff(end+1) = diff;
        else
            step_o2_diff(end+1) = diff;
        end
    end
end

% same for person x
step_x1_diff = [];
step_x2_diff = [];
real_x_idx = find(impacts(:,4) == 3 | impacts(:,4) == 4);
real_x_times = impacts(real_x_idx,1)./Fs_fsr;
for x = 2:length(real_x_times)
    diff = real_x_times(x)-real_x_times(x-1);
    if diff < 1.5 % not a turning point
        if impacts(real_x_idx(x),4) == 3 % left foot
            step_x1_diff(end+1) = diff;
        else
            step_x2_diff(end+1) = diff;
        end
    end
end

step_o1 = mean(step_o1_diff)
step_o2 = mean(step_o2_diff)
step_x1 = mean(step_x1_diff)
step_x2 = mean(step_x2_diff)

save(filename,'impactN','impacts','coordinates','acc_pk_idx','acc_pks',...
    'ID_labels','walk_edges','segments','step_o1','step_o2','step_x1','step_x2','-append')

%% sanity check walk edges

for i = 1:4
    figure; plot(fsrData(:,i))
    hold on
    impact_idx = find(impacts(:,4) == i);
    plot(impacts(impact_idx,1),fsrData(impacts(impact_idx,1),i),'rx')
end
%% wiener filter noise 5/17/22

% make sure these noise segments are pure noise
for s = 1:6
    figure;
    plot(datas(1:(Fs_pcb*15),s))
end

%% wiener filter noise cont. 5/17/22

wien_pcbD = zeros(length(pcbData),sensorN);
for s = 1:6
    noise_clip = datas(Fs_pcb*12:Fs_pcb*15,s);
    data = pcbData(:,s);
    [wiener_pcb,~] = WienerNoiseReduction([noise_clip;data],Fs_pcb,length(noise_clip));
    wiener_pcb = wiener_pcb(length(noise_clip)+1:end);
    wien_pcbD(:,s) = wiener_pcb;
end

% sanity check
figure; plot(pcbData(:,1))
title('Original')
figure; plot(wien_pcbD(:,1))
title('Wiener')

%% used to detect footfalls using cwt

freq_lower = 100; % Hz; frequency limits for cwt
freq_higher = 350;

final_estimates = [0,0,0]; % [impact time, segment ID num, likelihood score]
impact_thresh = 0.04; % (s) of how close an impact should be to the last

arrival_window = 1500; % 5000/Fs_pcb width of window to find arrival time
impact_time_thresh = 0.05; % if two impacts are 0.05s within each other, combine
correct_thresh = 0.05; % if predicted is within this thresh to real impact time, it's correct

success_rate = zeros(1,length(segments)-1);
failure_rate = zeros(1,length(segments)-1);
wrong_detect_rate = zeros(1,length(segments)-1);
big_failure_rate = zeros(1,length(segments)-1);

pkprom = [0.0015,0.0015,0.0015,0.0015,0.0015,0.0015];
pkheight = [0.006,0.016,0.0014,0.004,0.0127,0.009];

for w = 2:length(segments) % for each walking segment
    w
% for w=3:3
    start_time = impacts(segments(w-1),1)/Fs_fsr-0.15;
    stop_time = impacts(segments(w)-1,1)/Fs_fsr+0.15;
    
    start_idx_pcb = findTindex(start_time,pcbTime);
    stop_idx_pcb = findTindex(stop_time,pcbTime);
    
    impacttimes = impacts(segments(w-1):segments(w)-1,1)/Fs_fsr;
    
    estimated_impacts = [0,0]; % [impact time, sensor N], for each segment
    
    for s = 1:sensorN % for each sensor
        pcbclip = wien_pcbD(start_idx_pcb:stop_idx_pcb,s);
        [wt,f] = cwt(pcbclip,Fs_pcb); % uses default Morse wavelet
        valid_f_idx = find(f < freq_higher & f > freq_lower);
        cwt_freq = f(valid_f_idx);
        cwt_mag = abs(wt(valid_f_idx,:));
        sum_cwt = sum(cwt_mag,1);
        sum_smooth_cwt = movmean(sum_cwt, 1000);
        [pks, locs, ~,~] = findpeaks(sum_smooth_cwt,'MinPeakProminence',pkprom(s));%,'MinPeakHeight',pkheight(s));
        figure;
        findpeaks(sum_smooth_cwt,'MinPeakProminence',pkprom(s));%,'MinPeakHeight',pkheight(s))
        hold on
        impacttimes = impacts(segments(w-1):segments(w)-1,1)/Fs_fsr;
        plot((impacttimes-start_time)*Fs_pcb,0,'rx','MarkerSize',8)
        % new code 6/5/22
        for i = 1:length(locs)
            if i == 1
                starti = max(locs(i) - arrival_window, 1);
            else
                if locs(i) - arrival_window > locs(i-1) % isolated impact
                    starti = locs(i) - arrival_window;
                else % prior impact happens close before
                    [minval,minidx] = min(sum_smooth_cwt(locs(i-1):locs(i)));
                    starti = minidx + locs(i-1) - 1;
                end
            end
            window = sum_smooth_cwt(starti:locs(i));
            arrive_idx = aic_pick(window, 'to_peak')+starti - 1;
            add_array = [arrive_idx, s];
            estimated_impacts = [estimated_impacts; add_array];
            plot(arrive_idx,0,'bx') 
        end
        % old code, used until end of April 3 
%         for i = 1:length(locs)
%             if i == 1
%                 starti = max(locs(i) - arrival_window, 1);
%             else
%                 [minval, minidx] = min(sum_smooth_cwt(locs(i-1):locs(i)));
%                 starti = max(locs(i) - arrival_window, locs(i-1)+minidx - arrival_window/5);
%             end
%             if i ~= 1 & minval > 0.001
%                 arrive_idx = minidx + locs(i-1);
%             else
%                 window = sum_smooth_cwt(starti:locs(i));
%                 arrive_idx = aic_pick(window, 'to_peak')+starti;
%             end
%             add_array = [arrive_idx, s];
%             estimated_impacts = [estimated_impacts; add_array];
%             plot(arrive_idx,0,'bx') 
%         end
    end
    estimated_impacts(1,:) = []; % from initialization
    
    % combine estimated impacts so there's no repeat
    % new code: takes average of similar groups of impacts instead 6/5/22
    estimates = sort(estimated_impacts(:,1));
    current_impacts = [];
    scores = [];
    groupC = zeros(1,length(estimates));
    currval = estimates(1);
    groupC(1) = 1;
    thresh = 0.05*Fs_pcb;
    currCount = 1;
    for i = 2:length(estimates)
        if estimates(i) < currval + thresh
            groupC(i) = currCount;
        else
            currval = estimates(i);
            currCount = currCount + 1;
            groupC(i) = currCount;
        end
    end
    allcounts = unique(groupC);
    for i = 1:length(allcounts)
        groupidx = find(groupC == allcounts(i));
        current_impacts(end+1) = mean(estimates(groupidx));
        scores(end+1) = length(groupidx);
    end

    delete_pred = [];
    for i = 2:length(current_impacts)
        if current_impacts(i) - current_impacts(i-1) < thresh
            current_impacts(i) = mean(current_impacts(i-1:i));
            delete_pred(end+1) = i-1;
            scores(i) = scores(i) + scores(i-1);
        end
    end
    current_impacts(delete_pred) = [];
    scores(delete_pred) = [];
    
    % old code: takes sensor 1 for granted
%     sensor1_impact_idx = find(estimated_impacts(:,2) == 1);
%     % automatically add sensor1 impacts
%     current_impacts = estimated_impacts(sensor1_impact_idx,1);
%     scores = ones(1,length(current_impacts));
%     for s = 2:sensorN
%         sensor_impact_idx = find(estimated_impacts(:,2) == s);
%         for i = 1:length(sensor_impact_idx)
%             current_val = estimated_impacts(sensor_impact_idx(i),1);
%             current_difference = abs(current_impacts - current_val)/Fs_pcb;
%             near_value_idx = find(current_difference < impact_time_thresh);
%             % if current value is far enough away from recorded values
%             if isempty(near_value_idx)
%                 current_impacts(end+1) = current_val;
%                 scores(end+1) = 1;
%             else % repeat of an already recorded value
%                 scores(near_value_idx(1)) = scores(near_value_idx(1)) + 1;
%             end
%         end
%     end
    current_impacts = current_impacts + start_idx_pcb;
    elim = find(current_impacts < (impacttimes(1)-0.25)*Fs_pcb);
    current_impacts(elim) = [];
    scores(elim) = [];
    elim = find(current_impacts > (impacttimes(end)+0.25)*Fs_pcb);
    current_impacts(elim) = [];
    scores(elim) = [];
    % add current_impacts to final_estimates
    for i = 1:length(current_impacts)
        add_array = [current_impacts(i),w,scores(i)];
        final_estimates = [final_estimates; add_array];
    end
    
    figure;
    plot(impacttimes,0,'rx')
    hold on
    plot(current_impacts/Fs_pcb,1,'bx')
    ylim([-1 2])
    
    % calculating success rate
    current_impact_times = current_impacts./Fs_pcb;
    successN = length(impacttimes);
    success_count = 0;
    for i = 1:successN
        diff = abs(current_impact_times - impacttimes(i));
        if ~isempty(find(diff < correct_thresh))
            success_count = success_count + 1;
        end
    end
    estimatedN = length(current_impact_times);
    wrong_detect = 0;
    for i = 1:estimatedN
        diff = abs(impacttimes - current_impact_times(i));
        if isempty(find(diff < correct_thresh))
            wrong_detect = wrong_detect + 1;
        end
    end
    w
    success_rate(w-1) = success_count / successN
    failure_rate(w-1) = (successN - success_count)/successN
    wrong_detect_rate(w-1) = wrong_detect / estimatedN
    
    current_impact_times = current_impacts./Fs_pcb;
    successN = length(impacttimes);
    big_failure_count = 0;
    for i = 1:successN
        diff = abs(current_impact_times - impacttimes(i));
        if ~isempty(find(diff < 0.15))
            big_failure_count = big_failure_count + 1;
        end
    end
    big_failure_rate(w-1) = (successN - big_failure_count) / successN
end

final_estimates(1,:) = []; % from initialization
% UNCOMMENT to save variables
save(filename,'wien_pcbD','final_estimates','success_rate','failure_rate','wrong_detect_rate','big_failure_rate','-append')


%% using above, look at all combinations, find one with min uncertainty
% with LIMPS (same code as above section)
% instead of estimateID, use estimateID_LR

close all
final_estimates_labels = [0,0,0];
plotbool = true;
for s = 2:length(segments)
% for s = 9:length(segments)
% for s = 9:9
    s
    start_index = segments(s-1);
    stop_index = segments(s)-1;

    step_times_idx = find(final_estimates(:,2) == s);
    orig_step_times = final_estimates(step_times_idx,1)./Fs_pcb;
    step_times_scores = final_estimates(step_times_idx,3);
    seg_matrix = [orig_step_times,step_times_scores];
    sort_seg_matrix = sortrows(seg_matrix,1);
    step_times = sort_seg_matrix(:,1);
    step_times_scores = sort_seg_matrix(:,2);
    
    curr_ID_labels = ID_labels(start_index:stop_index);
    real_impact_times = impacts(start_index:stop_index,1)./Fs_fsr;
    
    if length(step_times) > 22 % uniqueperms can only handle 22 length
        diff = length(step_times) - 22;
        cut_start = floor(diff/2);
        cut_end = ceil(diff/2);
        step_times = step_times(cut_start+1:length(orig_step_times)-cut_end);
        step_times_scores = step_times_scores(cut_start+1:length(step_times_scores)-cut_end);
    end
    
    % trying to take out the recursion in the function
    [diff_array, label_array, estimateID] = notrecursive_stepID_limp_whole(curr_ID_labels,real_impact_times,step_times,step_times_scores,step_o1,step_o2,step_x1,step_x2,plotbool);
    est_o1 = find(label_array == 1);
    est_o2 = find(label_array == 2);
    est_x1 = find(label_array == 3);
    est_x2 = find(label_array == 4);

    bigidxone1 = find(diff_array(est_o1) > step_o1*1.58);
    bigidxone2 = find(diff_array(est_o2) > step_o2*1.58);
    bigidxtwo1 = find(diff_array(est_x1) > step_x1*1.58);
    bigidxtwo2 = find(diff_array(est_x2) > step_x2*1.58);
    bigidx = sort([est_o1(bigidxone1);est_o2(bigidxone2);est_x1(bigidxtwo1);est_x2(bigidxtwo2)]);
    smallidxone1 = find(diff_array(est_o1(2:end)) < step_o1/1.58);
    smallidxone2 = find(diff_array(est_o2(2:end)) < step_o2/1.58);
    smallidxtwo1 = find(diff_array(est_x1(2:end)) < step_x1/1.58);
    smallidxtwo2 = find(diff_array(est_x2(2:end)) < step_x2/1.58);
    smallidx = sort([est_o1(smallidxone1+1),est_o2(smallidxone2+1),est_x1(smallidxtwo1+1),est_x2(smallidxtwo2+1)]);
    whilecount = 0;
    while ~isempty(bigidx) | ~isempty(smallidx)
        whilecount = whilecount + 1
        if ~isempty(bigidx)
            [bad_idx,~] = max(bigidx);
            bad_idx % uncomment to debug
            curr_idx = bigidx(1);

            if label_array(curr_idx) == 1
                [~,change_idx] = min(abs(step_times(curr_idx)-step_times(1:curr_idx-1) - step_o1));
            elseif label_array(curr_idx) == 2
                [~,change_idx] = min(abs(step_times(curr_idx)-step_times(1:curr_idx-1) - step_o2));
            elseif label_array(curr_idx) == 3
                [~,change_idx] = min(abs(step_times(curr_idx)-step_times(1:curr_idx-1) - step_x1));
            else
                [~,change_idx] = min(abs(step_times(curr_idx)-step_times(1:curr_idx-1) - step_x2));
            end
            change_idx
            step_times = [step_times(1:change_idx);step_times(change_idx:end)];
            step_times_scores = [step_times_scores(1:change_idx);step_times_scores(change_idx:end)];
            [diff_array, label_array, estimateID] = notrecursive_stepID_limp_whole(curr_ID_labels,real_impact_times,step_times,step_times_scores,step_o1,step_o2,step_x1,step_x2,plotbool);
        else
            % delete element at smallidx or smallidx-1 depending on score
            curr_score = step_times_scores(smallidx(1));
            prev_score = step_times_scores(smallidx(1)-1);
            if curr_score > prev_score
                step_times = [step_times(1:smallidx(1)-2); step_times(smallidx(1):end)];
                step_times_scores(smallidx(1)-1) = [];
                smallidx(1)-1
    %             scores = [scores(1:smallidx(1)-2);scores(smallidx(1):end)];
            else
                step_times = [step_times(1:smallidx(1)-1); step_times(smallidx(1)+1:end)];
                step_times_scores(smallidx(1)) = [];
                smallidx(1)
    %             scores = [scores(1:smallidx(1)-1);scores(smallidx(1)+1:end)];
            end
            [diff_array, label_array, estimateID] = notrecursive_stepID_limp_whole(curr_ID_labels,real_impact_times,step_times,step_times_scores,step_o1,step_o2,step_x1,step_x2,plotbool);
        end
        est_o1 = find(label_array == 1);
        est_o2 = find(label_array == 2);
        est_x1 = find(label_array == 3);
        est_x2 = find(label_array == 4);

        bigidxone1 = find(diff_array(est_o1) > step_o1*1.58);
        bigidxone2 = find(diff_array(est_o2) > step_o2*1.58);
        bigidxtwo1 = find(diff_array(est_x1) > step_x1*1.58);
        bigidxtwo2 = find(diff_array(est_x2) > step_x2*1.58);
        bigidx = sort([est_o1(bigidxone1);est_o2(bigidxone2);est_x1(bigidxtwo1);est_x2(bigidxtwo2)]);
        smallidxone1 = find(diff_array(est_o1(2:end)) < step_o1/1.58);
        smallidxone2 = find(diff_array(est_o2(2:end)) < step_o2/1.58);
        smallidxtwo1 = find(diff_array(est_x1(2:end)) < step_x1/1.58);
        smallidxtwo2 = find(diff_array(est_x2(2:end)) < step_x2/1.58);
        smallidx = sort([est_o1(smallidxone1+1),est_o2(smallidxone2+1),est_x1(smallidxtwo1+1),est_x2(smallidxtwo2+1)]);
    end
    
    
    
    % tried with recursive_stepID_limp_whole.m, didn't work
%     if length(orig_step_times) > 20 % uniqueperms can only handle 22 length
%         diff = length(orig_step_times) - 20;
%         cut_start = floor(diff/2);
%         cut_end = ceil(diff/2);
%         orig_step_times = orig_step_times(cut_start+1:length(orig_step_times)-cut_end);
%         step_times_scores = step_times_scores(cut_start+1:length(step_times_scores)-cut_end);
%     end
%     
%     real_labels = ID_labels(start_index:stop_index);
%     real_impact_times = impacts(start_index:stop_index,1)./Fs_fsr;
%     
%     [estimateID_LR, step_times, scores] = recursive_stepID_limp_whole(real_labels, real_impact_times, orig_step_times,step_times_scores,step_o1,step_o2,step_x1,step_x2);
    
    
    
    
    
%     
% %     shorten the search by cutting segments in half
%     half_idx = floor(length(real_impact_times)/2);
%     if (real_impact_times(half_idx+1) - real_impact_times(half_idx)) < 0.1
%         half_idx = half_idx + 1;
%     end
%     half_time = real_impact_times(half_idx);
%     
%     % first half
%     real_idx = find(real_impact_times <= half_time);
%     step_times_idx = find(orig_step_times <= (half_time)); 
%     [estimateID_LR1, step_times1, scores1] = recursive_stepID_limp(real_labels(real_idx), real_impact_times(real_idx), orig_step_times(step_times_idx),step_times_scores(step_times_idx),step_o1,step_o2,step_x1,step_x2);
%     
%     title(sprintf('Final estimated ID labels %d segment'),s)
%     
%     % second half
%     real_idx = find(real_impact_times > half_time);
%     step_times_idx = find(orig_step_times > (half_time)); 
%     [estimateID_LR2, step_times2, scores2] = recursive_stepID_limp(real_labels(real_idx), real_impact_times(real_idx), orig_step_times(step_times_idx),step_times_scores(step_times_idx),step_o1,step_o2,step_x1,step_x2);
%     
%     title(sprintf('Final estimated ID labels %d segment'),s)
%     
%     % eliminate overlap between the two halves
%     delete_idx2 = [];
%     for i = 1:length(step_times2)
%         curr_step_time = step_times2(i);
%         same_idx = find(step_times1 == curr_step_time);
%         if ~isempty(same_idx)
%             for j = 1:length(same_idx)
%                 curr_id = estimateID_LR2(i);
%                 same_id = estimateID_LR1(same_idx(j));
%                 if curr_id == 1 | curr_id == 2
%                     if same_id == 1 | same_id == 2
%                         % overlapping impact
%                         delete_idx2(end+1) = i;
%                     end
%                 else
%                     if same_id == 3 | same_id == 4
%                         % overlapping impact
%                         delete_idx2(end+1) = i;
%                     end
%                 end
%             end
%         end
%     end
%     step_times2(delete_idx2) = [];
%     estimateID_LR2(delete_idx2) = [];
%     scores2(delete_idx2)= [];
%     
%     % combine the two halves, add last edits
%     [estimateID, step_times] = recursive_stepID_noperms_limp(real_labels,real_impact_times,[estimateID_LR1,estimateID_LR2],[step_times1;step_times2],[scores1;scores2],step_o1,step_o2,step_x1,step_x2);
%     
    for i = 1:length(step_times)
        add_array = [step_times(i),label_array(i),s];
        final_estimates_labels = [final_estimates_labels; add_array];
    end
    
end

final_estimates_labels(1,:) = []; % from initialization

figure;
real_o = find(ID_labels == 1);
real_x = find(ID_labels == 2);
plot(impacts(real_o,1)./Fs_fsr,0,'bo')
hold on
plot(impacts(real_x,1)./Fs_fsr,0,'bx')
ylim([-1 1])

est_o = find(final_estimates_labels(:,2) == 1 | final_estimates_labels(:,2) == 2);
est_x = find(final_estimates_labels(:,2) == 3 | final_estimates_labels(:,2) == 4);

plot(final_estimates_labels(est_o,1),0.5,'ro')
plot(final_estimates_labels(est_x,1),0.5,'rx')

save(filename,'final_estimates_labels','-append')

%% remove start/end of real segments to match final estimates 5/21/22

clean_impacts = impacts;
clean_ID_labels = ID_labels;
seg_val = unique(final_estimates_labels(:,3));
walk_edges_starts = find(walk_edges == -1);
walk_edges_ends = find(walk_edges == 1);
delete_idx = [];
for i = 1:length(seg_val)
    seg_idx = find(final_estimates_labels(:,3) == seg_val(i));
    first_val = final_estimates_labels(seg_idx(1),1);
    last_val = final_estimates_labels(seg_idx(end),1);
    
    impact_inseg = impacts(walk_edges_starts(seg_val(i)-1):walk_edges_ends(seg_val(i)-1),1)./Fs_fsr;
    curr_delete_idx = find(impact_inseg < (first_val - 0.15) | impact_inseg > (last_val + 0.15));
    delete_idx = [delete_idx;curr_delete_idx + walk_edges_starts(seg_val(i)-1) - 1];
end

clean_impacts(delete_idx,:) = [];
clean_ID_labels(delete_idx) = [];

figure;
real_o = find(clean_ID_labels == 1);
real_x = find(clean_ID_labels == 2);
plot(clean_impacts(real_o,1)./Fs_fsr,0,'bo')
hold on
plot(clean_impacts(real_x,1)./Fs_fsr,0,'bx')
ylim([-1 1])

est_o = find(final_estimates_labels(:,2) == 1 | final_estimates_labels(:,2) == 2);
est_x = find(final_estimates_labels(:,2) == 3 | final_estimates_labels(:,2) == 4);

plot(final_estimates_labels(est_o,1),0.5,'ro')
plot(final_estimates_labels(est_x,1),0.5,'rx')

figure; plot(clean_impacts(real_o,1)./Fs_fsr,real_o,'bo')
hold on
plot(clean_impacts(real_x,1)./Fs_fsr,real_x,'bx')

figure; plot(final_estimates_labels(est_o,1),est_o,'ro')
hold on
plot(final_estimates_labels(est_x,1),est_x,'rx')

%% manually remove beg/end of segments that are extra impacts (4/25/22) 
% clean_... means start/end of segments that don't match have been removed

% same as above but with index
remove_est_idx = [109,110,111,129,192]; % seconds
remove_real_idx = [18,40,42,44,101,121,122,123,125,127,145,209]; % seconds

clean_final_estimates_labels = final_estimates_labels;
clean_final_estimates_labels(remove_est_idx,:) = [];


clean_impacts(remove_real_idx,:) = [];
clean_ID_labels(remove_real_idx) = [];

% plot results
figure;
real_o = find(clean_ID_labels == 1);
real_x = find(clean_ID_labels == 2);
plot(clean_impacts(real_o,1)./Fs_fsr,0,'bo')
hold on
plot(clean_impacts(real_x,1)./Fs_fsr,0,'bx')
ylim([-1 1])

est_o = find(clean_final_estimates_labels(:,2) == 1 | clean_final_estimates_labels(:,2) == 2);
est_x = find(clean_final_estimates_labels(:,2) == 3 | clean_final_estimates_labels(:,2) == 4);

plot(clean_final_estimates_labels(est_o,1),0.5,'ro')
plot(clean_final_estimates_labels(est_x,1),0.5,'rx')

save(filename,'clean_final_estimates_labels','clean_impacts','clean_ID_labels','-append')


%% add overlapping impacts when alg didn't pick up bc start/end of segment 5/21/22

add_impacts_beg = [1.8497,23.5991,62.9905,70.6383]; % (s); at beg of episode; these need to be in proper order - later ones first
add_impacts_end = [3.83492,4.31742,11.9941,27.7414,43.3963,67.1973,74.7644,83.103]; % (s); at end of episode

for i = 1:length(add_impacts_beg)
    [~,idx] = min(abs(clean_final_estimates_labels(:,1) - add_impacts_beg(i)));
    curr_label = clean_final_estimates_labels(idx,2);
    new_label = 0;
    if curr_label == 1 | curr_label == 2
        next_idx = idx+1;
        next_label = clean_final_estimates_labels(next_idx,2);
        while(next_label ~= 3 & next_label ~=4)
            next_idx = next_idx+1;
            next_label = clean_final_estimates_labels(next_idx,2);
        end
        if next_label == 3
            new_label = 4;
        else
            new_label = 3;
        end
    else
        next_idx = idx+1;
        next_label = clean_final_estimates_labels(next_idx,2);
        while(next_label ~= 1 & next_label ~=2)
            next_idx = next_idx+1;
            next_label = clean_final_estimates_labels(next_idx,2);
        end
        if next_label == 1
            new_label = 2;
        else
            new_label = 1;
        end
    end
    insert_arr = [clean_final_estimates_labels(idx,1),new_label,clean_final_estimates_labels(idx,3)];
    clean_final_estimates_labels = [clean_final_estimates_labels(1:idx,:);insert_arr;clean_final_estimates_labels(idx+1:end,:)];
end

for i = 1:length(add_impacts_end)
    [~,idx] = min(abs(clean_final_estimates_labels(:,1) - add_impacts_end(i)));
    curr_label = clean_final_estimates_labels(idx,2);
    new_label = 0;
    prev_idx = idx-1;
    prev_label = clean_final_estimates_labels(prev_idx,2);
    if curr_label == 1 | curr_label == 2
        while(prev_label ~= 3 & prev_label ~=4)
            prev_idx = prev_idx-1;
            prev_label = clean_final_estimates_labels(prev_idx,2);
        end
        if prev_label == 3
            new_label = 4;
        else
            new_label = 3;
        end
    else
        while(prev_label ~= 1 & prev_label ~=2)
            prev_idx = prev_idx-1;
            prev_label = clean_final_estimates_labels(prev_idx,2);
        end
        if prev_label == 1
            new_label = 2;
        else
            new_label = 1;
        end
    end
    insert_arr = [clean_final_estimates_labels(idx,1),new_label,clean_final_estimates_labels(idx,3)];
    clean_final_estimates_labels = [clean_final_estimates_labels(1:idx,:);insert_arr;clean_final_estimates_labels(idx+1:end,:)];
end

figure;
real_o = find(clean_ID_labels == 1);
real_x = find(clean_ID_labels == 2);
plot(clean_impacts(real_o,1)./Fs_fsr,0,'bo')
hold on
plot(clean_impacts(real_x,1)./Fs_fsr,0,'bx')
ylim([-1 1])

est_o = find(clean_final_estimates_labels(:,2) == 1 | clean_final_estimates_labels(:,2) == 2);
est_x = find(clean_final_estimates_labels(:,2) == 3 | clean_final_estimates_labels(:,2) == 4);

plot(clean_final_estimates_labels(est_o,1),0.5,'ro')
plot(clean_final_estimates_labels(est_x,1),0.5,'rx')

save(filename,'clean_final_estimates_labels','-append')

%% evaluate clean_final_estimates_labels (4/22/22)

segmentN = length(segments)-1;
false_pos = zeros(1,segmentN);
false_neg = zeros(1,segmentN);
true_pos = zeros(1,segmentN);
rmse = zeros(1,segmentN);
error_thresh = 0.05;

big_false_pos = zeros(1,segmentN);
big_false_neg = zeros(1,segmentN);
big_true_pos = zeros(1,segmentN);
big_rmse = zeros(1,segmentN);
big_error_thresh = 0.1;

for s = 2:length(segments)
    
    curr_false_pos = 0;
    curr_false_neg = 0;
    curr_true_pos = 0;
    curr_rmse = 0;
    
    big_curr_false_pos = 0;
    big_curr_false_neg = 0;
    big_curr_true_pos = 0;
    big_curr_rmse = 0;
    
    start_index = segments(s-1);
    stop_index = segments(s)-1;
    start_time = impacts(start_index,1)./Fs_fsr;
    stop_time = impacts(stop_index,1)./Fs_fsr;
    
    clean_idx = find(clean_impacts(:,1)./Fs_fsr >= start_time & clean_impacts(:,1)./Fs_fsr <= stop_time);
    
    real_labels = clean_ID_labels(clean_idx);
    real_times = clean_impacts(clean_idx,1)./Fs_fsr;
    
    final_estimates_idx = find(clean_final_estimates_labels(:,3) == s); % find labels belonging to this segment number
    estimate_times = clean_final_estimates_labels(final_estimates_idx,1);
    estimate_labels = clean_final_estimates_labels(final_estimates_idx,2);
    % correct these for 1,2 = 2, 3,4 = 2
    for i = 1:length(estimate_labels)
        curr_label = estimate_labels(i);
        if curr_label == 2
            estimate_labels(i) = 1;
        elseif curr_label == 3 | curr_label == 4
            estimate_labels(i) = 2;
        end
    end
    
    realN = length(real_times);
    estN = length(estimate_times);
    
    % compare each real time with estimated times to find true+ & false-
    for t = 1:realN
        curr_time = real_times(t); % real impact time
        curr_label = real_labels(t);
        diff = abs(estimate_times - curr_time);
        diff_idx = find(diff < error_thresh);
        if ~isempty(diff_idx) % if a estimated impact exists close to the real impact time
            est_idx = find(estimate_labels(diff_idx)==curr_label);
            if ~isempty(est_idx) % if the close impact is labeled correctly
                curr_true_pos = curr_true_pos + 1;
                est_times = estimate_times(diff_idx);
                est_time = est_times(est_idx);
                curr_rmse = curr_rmse + (curr_time - est_time)^2; % rmse calculation
            else
                curr_false_neg = curr_false_neg + 1;
            end
        else
            curr_false_neg = curr_false_neg + 1;
        end
        
        diff_idx = find(diff < big_error_thresh);
        if ~isempty(diff_idx) % if a estimated impact exists close to the real impact time
            est_idx = find(estimate_labels(diff_idx)==curr_label);
            if ~isempty(est_idx) % if the close impact is labeled correctly
                [~,thisidx] = min(diff(est_idx));
                est_idx = est_idx(thisidx);
                big_curr_true_pos = big_curr_true_pos + 1;
                est_times = estimate_times(diff_idx);
                est_time = est_times(est_idx);
                big_curr_rmse = big_curr_rmse + (curr_time - est_time)^2; % rmse calculation
            else
                big_curr_false_neg = big_curr_false_neg + 1;
            end
        else
            big_curr_false_neg = big_curr_false_neg + 1;
        end
    end
    
    % compare each est time with real times to find false+
    for t = 1:estN
        curr_est_time = estimate_times(t); % estimated impact time
        curr_est_label = estimate_labels(t);
        diff = abs(real_times - curr_est_time);
        diff_idx = find(diff < error_thresh);
        if ~isempty(diff_idx)
            real_idx = find(real_labels(diff_idx)==curr_est_label);
            if isempty(real_idx)
                curr_false_pos = curr_false_pos + 1;
            end
        else
            curr_false_pos = curr_false_pos + 1;
        end
        
        diff_idx = find(diff < big_error_thresh);
        if ~isempty(diff_idx)
            real_idx = find(real_labels(diff_idx)==curr_est_label);
            if isempty(real_idx)
                big_curr_false_pos = big_curr_false_pos + 1;
            end
        else
            big_curr_false_pos = big_curr_false_pos + 1;
        end
    end
    
    % store
    false_pos(s-1) = curr_false_pos/estN;
    false_neg(s-1) = curr_false_neg/realN;
    true_pos(s-1) = curr_true_pos/realN;
    rmse(s-1) = sqrt(curr_rmse/curr_true_pos);
    
    big_false_pos(s-1) = big_curr_false_pos/estN;
    big_false_neg(s-1) = big_curr_false_neg/realN;
    big_true_pos(s-1) = big_curr_true_pos/realN;
    big_rmse(s-1) = sqrt(big_curr_rmse/big_curr_true_pos);
end

true_pos
save(filename,'false_pos','false_neg','true_pos','rmse','big_false_pos','big_false_neg','big_true_pos','big_rmse','-append')

%% delete segment 1 5/25/22

i = 4; % this is +1 from count of segment
est_idx = find(clean_final_estimates_labels(:,3) == i);
clean_final_estimates_labels(est_idx,:) = [];
start_idx = find(clean_final_estimates_labels(:,3) == i+1);
start_idx = start_idx(1);
last_idx = find(clean_final_estimates_labels(:,3) == i-1);
if isempty(last_idx)
    real_idx = find(clean_impacts(:,1)./Fs_fsr < clean_final_estimates_labels(start_idx,1) - 0.15 & clean_impacts(:,1)./Fs_fsr > 0);
else
    last_idx = last_idx(end);
    real_idx = find(clean_impacts(:,1)./Fs_fsr < clean_final_estimates_labels(start_idx,1) - 0.15 & clean_impacts(:,1)./Fs_fsr > clean_final_estimates_labels(last_idx,1) + 0.15);
end
clean_impacts(real_idx,:) = [];
clean_ID_labels(real_idx) = [];

save(filename,'clean_final_estimates_labels','clean_impacts','clean_ID_labels','-append')
% run above section again to evaluate

%% perform GMM on these estimated step times for one intervention 5/19/22

estimated_scaled_means = zeros(2,2);
scaled_means = zeros(2,2);
real_means = zeros(2,2);
real_stds = zeros(2,2);
for p = 1:2 % for each person
    idx = find(floor(clean_final_estimates_labels(:,2)./3) + 1 == p);
    estimated_impact_times = clean_final_estimates_labels(idx,1);

    % calculate all step times from impact times
    step_times = [];
    for x = 2:length(estimated_impact_times)
        curr_diff = estimated_impact_times(x) - estimated_impact_times(x-1);
        if curr_diff < 1.5 % not a turning point
            step_times(end+1) = curr_diff;
        end
    end

    % GM analysis from estimated step times
    GM = fitgmdist(transpose(step_times),2,'RegularizationValue',0.000001);
    proportion = GM.ComponentProportion;
    mu = GM.mu;
    [bigprop,~] = max(proportion);
    person_idx = (p-1)*2+1;
    if (1-bigprop) < abs(0.5-bigprop)
        mean1 = mu(1) + (mu(2)-mu(1))*(1- (proportion(1))^2);
        mean2 = mu(2) + (mu(1)-mu(2))*(1- (proportion(2))^2);
    else
        mean1 = mu(1) + (mu(2)-mu(1))*abs(0.5-proportion(1));
        mean2 = mu(2) + (mu(1)-mu(2))*abs(0.5-proportion(2));
    end
    estimated_scaled_means(p,:) = [min(mean1,mean2), max(mean1,mean2)];
        
    % calculate real step times
    left_right_diff = [];
    right_left_diff = [];
    real_differences = [];
    real_idx = find(clean_impacts(:,4) == person_idx | clean_impacts(:,4) == (person_idx+1));
    real_impact_times = clean_impacts(real_idx,1)./Fs_fsr;
    for x = 2:length(real_impact_times)
        real_diff = real_impact_times(x)-real_impact_times(x-1);
        if real_diff < 1.5 % not a turning point
            real_differences(end+1) = real_diff;
            if mod(clean_impacts(real_idx(x),4),2) == 0 % right foot
                left_right_diff(end+1) = real_diff;
            else
                right_left_diff(end+1) = real_diff;
            end
        end
    end

    % GM analysis from real step times
    GM = fitgmdist(transpose(real_differences),2,'RegularizationValue',0.000001);
    proportion = GM.ComponentProportion;
    mu = GM.mu;
    [bigprop,~] = max(proportion);
    if (1-bigprop) < abs(0.5-bigprop)
        mean1 = mu(1) + (mu(2)-mu(1))*(1- (proportion(1))^2);
        mean2 = mu(2) + (mu(1)-mu(2))*(1- (proportion(2))^2);
    else
        mean1 = mu(1) + (mu(2)-mu(1))*abs(0.5-proportion(1));
        mean2 = mu(2) + (mu(1)-mu(2))*abs(0.5-proportion(2));
    end
    scaled_means(p,:) = [min(mean1,mean2),max(mean1,mean2)];

    % real step times average & std
    mean1 = mean(left_right_diff);
    mean2 = mean(right_left_diff);
    real_means(p,:) = [min(mean1,mean2),max(mean1,mean2)];
    real_stds(p,:) = [std(left_right_diff),std(right_left_diff)];
end

estimated_scaled_means
scaled_means
real_means
save(filename,'estimated_scaled_means','scaled_means','real_means','real_stds','-append')

%% perform GMM on these estimated step times (4/26/22)
% 
% estimated_scaled_means = zeros(6,4);
% scaled_means = zeros(6,4);
% real_means = zeros(6,4);
% real_std = zeros(6,4);
% 
% intervention_arr = ['regular1','limp1','limp2','weight1','weight2','regular2'];
% file_root = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Jenny 1\ProcessedData\both_';
% 
% for i = 1:6
%     load([file_root,intervention_arr(i)])
%     for p = 1:2 % for each person
%         idx = find(clean_final_estimates_labels(:,2) == p);
%         estimated_impact_times = clean_final_estimates_labels(idx,1);
%         
%         % calculate all step times from impact times
%         step_times = [];
%         for x = 2:length(estimated_impact_times)
%             curr_diff = estimated_impact_times(x) - estimated_impact_times(x-1);
%             if curr_diff < 1.5 % not a turning point
%                 step_times(end+1) = curr_diff;
%             end
%         end
%         
%         % GM analysis from estimated step times
%         GM = fitgmdist(transpose(step_times),2,'RegularizationValue',0.000001);
%         proportion = GM.ComponentProportion;
%         mu = GM.mu;
%         [bigprop,~] = max(proportion);
%         person_idx = (p-1)*2+1;
%         if (1-bigprop) < abs(0.5-bigprop)
%             mean1 = mu(1) + (mu(2)-mu(1))*(1- (proportion(1))^2);
%             mean2 = mu(2) + (mu(1)-mu(2))*(1- (proportion(2))^2);
%         else
%             mean1 = mu(1) + (mu(2)-mu(1))*abs(0.5-proportion(1));
%             mean2 = mu(2) + (mu(1)-mu(2))*abs(0.5-proportion(2));
%         end
%         estimated_scaled_means(i,person_idx) = min(mean1,mean2);
%         estimated_scaled_means(i,person_idx+1) = max(mean1,mean2);
%         
%         % calculate real step times
%         left_right_diff = [];
%         right_left_diff = [];
%         real_differences = [];
%         real_idx = find(clean_impacts(:,4) == person_idx | clean_impacts(:,4) == (person_idx+1));
%         real_impact_times = clean_impacts(real_idx,1)./Fs_fsr;
%         for x = 2:length(real_impact_times)
%             real_diff = real_impact_times(x)-real_impact_times(x-1);
%             if real_diff < 1.5 % not a turning point
%                 real_differences(end+1) = real_diff;
%                 if mod(clean_impacts(real_idx(x),4),2) == 0 % right foot
%                     left_right_diff(end+1) = real_diff;
%                 else
%                     right_left_diff(end+1) = real_diff;
%                 end
%             end
%         end
%         
%         % GM analysis from real step times
%         GM = fitgmdist(transpose(real_differences),2,'RegularizationValue',0.000001);
%         proportion = GM.ComponentProportion;
%         mu = GM.mu;
%         [bigprop,~] = max(proportion);
%         if (1-bigprop) < abs(0.5-bigprop)
%             mean1 = mu(1) + (mu(2)-mu(1))*(1- (proportion(1))^2);
%             mean2 = mu(2) + (mu(1)-mu(2))*(1- (proportion(2))^2);
%         else
%             mean1 = mu(1) + (mu(2)-mu(1))*abs(0.5-proportion(1));
%             mean2 = mu(2) + (mu(1)-mu(2))*abs(0.5-proportion(2));
%         end
%         scaled_means(i,person_idx) = min(mean1,mean2);
%         scaled_means(i,person_idx+1) = max(mean1,mean2);
%         
%         % real step times average & std
%         mean1 = mean(left_right_diff);
%         mean2 = mean(right_left_diff);
%         real_means(i,person_idx) = min(mean1,mean2);
%         real_means(i,person_idx+1) = max(mean1,mean2);
%         real_std(i,person_idx) = std(left_right_diff);
%         real_std(i,person_idx+1) = std(right_left_diff);
%     end
% end
% 
% % save(filename,'estimated_scaled_means','scaled_means','real_means','real_std','-append')


%% using above, look at all combinations, find one with min uncertainty
% SKIP THIS!!!!! (doesn't consider limps)

% close all
% final_estimates_labels = [0,0,0];
% % for s = 2:length(segments)
% for s = 2:2
%     start_index = segments(s-1);
%     stop_index = segments(s)-1;
% 
%     step_times_idx = find(final_estimates(:,2) == s);
%     orig_step_times = final_estimates(step_times_idx,1)./Fs_pcb;
%     
%     real_labels = ID_labels(start_index:stop_index);
%     real_impact_times = impacts(start_index:stop_index,1)./Fs_fsr;
%     
% %     shorten the search by cutting segments in half
%     half_idx = floor(length(real_impact_times)/2);
%     if (real_impact_times(half_idx+1) - real_impact_times(half_idx)) < 0.1
%         half_idx = half_idx + 1;
%     end
%     half_time = real_impact_times(half_idx);
%     
%     real_idx = find(real_impact_times <= half_time);
%     step_times_idx = find(orig_step_times <= (half_time+0.1)); % includes a buffer for time cutoff
%     [estimateID1, step_times1] = recursive_stepID(real_labels(real_idx), real_impact_times(real_idx), sort(orig_step_times(step_times_idx)),step_o,step_x);
%     
%     title(sprintf('Final estimated ID labels %d segment'),s)
%     
%     % second half
%     real_idx = find(real_impact_times > half_time);
%     step_times_idx = find(orig_step_times > (half_time-0.1)); % includes a buffer for time cutoff
%     [estimateID2, step_times2] = recursive_stepID(real_labels(real_idx), real_impact_times(real_idx), sort(orig_step_times(step_times_idx)),step_o,step_x);
%     
%     title(sprintf('Final estimated ID labels %d segment'),s)
%     
%     % combine the two halves, add last edits
%     [estimateID, step_times] = recursive_stepID_noperms(real_labels, real_impact_times, [estimateID1,estimateID2], [step_times1;step_times2], step_o, step_x);
%     
%     for i = 1:length(step_times)
%         add_array = [step_times(i),estimateID(i),s];
%         final_estimates_labels = [final_estimates_labels; add_array];
%     end
%     
% end
% 
% final_estimates_labels(1,:) = []; % from initialization
% 
% figure;
% real_o = find(ID_labels == 1);
% real_x = find(ID_labels == 2);
% plot(impacts(real_o,1)./Fs_fsr,0,'bo')
% hold on
% plot(impacts(real_x,1)./Fs_fsr,0,'bx')
% ylim([-1 1])
% 
% est_o = find(final_estimates_labels(:,2) == 1);
% est_x = find(final_estimates_labels(:,2) == 2);
% 
% plot(final_estimates_labels(est_o,1),0.5,'ro')
% plot(final_estimates_labels(est_x,1),0.5,'rx')

% save(filename,'final_estimates_labels','-append')




































