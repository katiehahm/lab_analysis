%% 4/18/22 experiment 4 all processing, do it section by section
clear all
close all
filepath = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Uriel 2\';
% ################# change ########################
intervention = 'both_weight1'; % regular1 limp1 lim2 weight1 weight2 regular2
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

start_time = 17.3; % change ##
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


last_time = 172 - start_time; % change ##
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

%%
% clean and filter pcb data
filt_pcbD = lpf_data(pcbData);

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

%% fix small errors in impacts 11/5/21
heel_start_wrong = [8458]; % these need to be same length
heel_start_right = [8544];

heel_pk_wrong = []; % index, these need to be same length
heel_pk_right = [];

delete_pks = []; % index of peaks to delete

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
processedfilepath = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Jenny 1\ProcessedData\';
filename = [processedfilepath,intervention];

save(filename,'pcbTime','filt_pcbD','datas','pcbData', ...
    'fsrTime','fsrData','impacts','Fs_pcb','Fs_fsr','Fs_acc','Fs_trigger', ...
    'acc_pks','acc_pk_idx','accTime','accData',...
    'mocapT','allmocap','coordinates','sensorN')
disp(append("Saved as ", filename))

%% FOOTFALL DETECTION

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

% find average step time for each person
% # CHANGE FOR LIMPING ###################
ID_o = 0;
ID_o_count = 0;
last_o_idx = 0;
ID_x = 0;
ID_x_count = 0;
last_x_idx = 0;
for i = 1:2
    if ID_labels(i) == 1
        last_o_idx = i;
    else
        last_x_idx = i;
    end
end
for i = 3:impactN % assumes first two impacts are from diff person
    if walk_edges(i) == -1 % start of new walking segment
        % first step of segment
        if ID_labels(i) == 1
            last_o_idx = i;
        else
            last_x_idx = i;
        end
        % second step of segment
        if ID_labels(i+1) == 1
            last_o_idx = i+1;
        else
            last_x_idx = i+1;
        end
        i = i + 2; % advance to 3rd step in segment
    else
        if ID_labels(i) == 1
            difference = impacts(i,1) - impacts(last_o_idx,1);
            ID_o = ID_o + difference/Fs_fsr;
            ID_o_count = ID_o_count + 1;
            last_o_idx = i;
        else
            difference = impacts(i,1) - impacts(last_x_idx,1);
            ID_x = ID_x + difference/Fs_fsr;
            ID_x_count = ID_x_count + 1;
            last_x_idx = i;
        end
    end
end

step_o = ID_o/ID_o_count
step_x = ID_x/ID_x_count

save(filename,'impactN','impacts','coordinates','acc_pk_idx','acc_pks',...
    'ID_labels','walk_edges','segments','step_o','step_x','-append')

%% sanity check walk edges
figure; plot(fsrData(:,1))
hold on
impact_1_idx = find(impacts(:,4) == 1);
plot(impacts(impact_1_idx,1),fsrData(impacts(impact_1_idx,1),1),'rx')

%% used to detect footfalls using cwt

freq_lower = 100; % Hz; frequency limits for cwt
freq_higher = 350;

final_estimates = [0,0]; % [impact time, segment ID num]
impact_thresh = 0.04; % (s) of how close an impact should be to the last

arrival_window = 3000; % 5000/Fs_pcb width of window to find arrival time
impact_time_thresh = 0.1; % if two impacts are 0.1s within each other, combine
correct_thresh = 0.05; % if predicted is within this thresh to real impact time, it's correct

% for w = 2:length(segments) % for each walking segment
for w=16:16
    start_time = impacts(segments(w-1),1)/Fs_fsr-0.15;
    stop_time = impacts(segments(w)-1,1)/Fs_fsr+0.15;
    
    start_idx_pcb = findTindex(start_time,pcbTime);
    stop_idx_pcb = findTindex(stop_time,pcbTime);
    
    impacttimes = impacts(segments(w-1):segments(w)-1,1)/Fs_fsr;
    
    estimated_impacts = [0,0]; % [impact time, sensor N], for each segment
    
    for s = 1:sensorN % for each sensor
        pcbclip = pcbData(start_idx_pcb:stop_idx_pcb,s);
        [wt,f] = cwt(pcbclip,Fs_pcb); % uses default Morse wavelet
        valid_f_idx = find(f < freq_higher & f > freq_lower);
        cwt_freq = f(valid_f_idx);
        cwt_mag = abs(wt(valid_f_idx,:));
        sum_cwt = sum(cwt_mag,1);
        sum_smooth_cwt = movmean(sum_cwt, 800);
        [pks, locs, ~,~] = findpeaks(sum_smooth_cwt, 'MinPeakProminence',0.0015);
        figure;
        findpeaks(sum_smooth_cwt, 'MinPeakProminence',0.0015)
        hold on
        impacttimes = impacts(segments(w-1):segments(w)-1,1)/Fs_fsr;
        plot((impacttimes-start_time)*Fs_pcb,0,'rx','MarkerSize',8)
        for i = 1:length(locs)
            if i == 1
                starti = max(locs(i) - arrival_window, 1);
            else
                [minval, minidx] = min(sum_smooth_cwt(locs(i-1):locs(i)));
                starti = max(locs(i) - arrival_window, locs(i-1)+minidx - arrival_window/5);
            end
            if i ~= 1 & minval > 0.001
                arrive_idx = minidx + locs(i-1);
            else
                window = sum_smooth_cwt(starti:locs(i));
                arrive_idx = aic_pick(window, 'to_peak')+starti;
            end
            add_array = [arrive_idx, s];
            estimated_impacts = [estimated_impacts; add_array];
            plot(arrive_idx,0,'bx') 
        end
        
    end
    
    estimated_impacts(1,:) = []; % from initialization
    sensor1_impact_idx = find(estimated_impacts(:,2) == 1);
    % automatically add sensor1 impacts
    current_impacts = estimated_impacts(sensor1_impact_idx,1);
    for s = 2:sensorN
        sensor_impact_idx = find(estimated_impacts(:,2) == s);
        for i = 1:length(sensor_impact_idx)
            current_val = estimated_impacts(sensor_impact_idx(i),1);
            current_difference = abs(current_impacts - current_val)/Fs_pcb;
            % if current value is far enough away from recorded values
            if isempty(find(current_difference < impact_time_thresh))
                current_impacts(end+1) = current_val;
            end
        end
    end
    current_impacts = current_impacts + start_idx_pcb;
    elim = find(current_impacts < (impacttimes(1)-0.25)*Fs_pcb);
    current_impacts(elim) = [];
    elim = find(current_impacts > (impacttimes(end)+0.25)*Fs_pcb);
    current_impacts(elim) = [];
    % add current_impacts to final_estimates
    for i = 1:length(current_impacts)
        add_array = [current_impacts(i),w];
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
    success_rate = success_count / successN
    failure_rate = (successN - success_count)/successN
    wrong_detect_rate = wrong_detect / estimatedN
    
    current_impact_times = current_impacts./Fs_pcb;
    successN = length(impacttimes);
    big_failure_count = 0;
    for i = 1:successN
        diff = abs(current_impact_times - impacttimes(i));
        if ~isempty(find(diff < 0.15))
            big_failure_count = big_failure_count + 1;
        end
    end
    big_failure_rate = (successN - big_failure_count) / successN
end

final_estimates(1,:) = []; % from initialization
% save(filename,'final_estimates','-append')
% edited for this file up to here, lower code is still just copied
% ################################################
%% using above, look at all combinations, find one with min uncertainty
% will take a long time to run

% close all
final_estimates_labels = [0,0];
for s = 2:length(segments)
% for s = 11:11
    start_index = segments(s-1);
    stop_index = segments(s)-1;

    step_times_idx = find(final_estimates(:,2) == s);
    orig_step_times = final_estimates(step_times_idx,1)./Fs_pcb;
    
    curr_ID_labels = ID_labels(start_index:stop_index);
    real_impact_times = impacts(start_index:stop_index,1)./Fs_fsr;
    
%     shorten the search by cutting segments in half
    half_idx = floor(length(real_impact_times)/2);
    if (real_impact_times(half_idx+1) - real_impact_times(half_idx)) < 0.1
        half_idx = half_idx + 1;
    end
    half_time = real_impact_times(half_idx);
    
    real_idx = find(real_impact_times <= half_time);
    step_times_idx = find(orig_step_times <= (half_time+0.1)); % includes a buffer for time cutoff
    [estimateID1, step_times1] = recursive_stepID(curr_ID_labels(real_idx), real_impact_times(real_idx), sort(orig_step_times(step_times_idx)),step_o,step_x);

%     [estimateID, step_times] = recursive_stepID(curr_ID_labels, real_impact_times, sort(step_times),step_o,step_x);

    
    title(sprintf('Final estimated ID labels %d segment'),s)
    
    for i = 1:length(step_times)
        add_array = [step_times(i),estimateID(i)];
        final_estimates_labels = [final_estimates_labels; add_array];
    end
    
    % second half
    real_idx = find(real_impact_times > half_time);
    step_times_idx = find(orig_step_times > (half_time-0.1)); % includes a buffer for time cutoff
    [estimateID2, step_times2] = recursive_stepID(curr_ID_labels(real_idx), real_impact_times(real_idx), sort(orig_step_times(step_times_idx)),step_o,step_x);

%     [estimateID, step_times] = recursive_stepID(curr_ID_labels, real_impact_times, sort(step_times),step_o,step_x);

    
    title(sprintf('Final estimated ID labels %d segment'),s)
    
    [estimateID, step_times] = recursive_stepID_noperms(curr_ID_labels, real_impact_times, [estimateID1,estimateID2], [step_times1;step_times2], step_o, step_x);
    
    for i = 1:length(step_times)
        add_array = [step_times(i),estimateID(i)];
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

est_o = find(final_estimates_labels(:,2) == 1);
est_x = find(final_estimates_labels(:,2) == 2);

plot(final_estimates_labels(est_o,1),0.5,'ro')
plot(final_estimates_labels(est_x,1),0.5,'rx')

% save(filename,'final_estimates_labels','-append')

%% evaluate final_estimates_labels (4/22/22)

segmentN = length(segments)-1;
false_pos = zeros(1,segmentN);
false_neg = zeros(1,segmentN);
true_pos = zeros(1,segmentN);
rmse = zeros(1,segmentN);

for i=2:length(segments)
    
    curr_false_pos = 0;
    curr_false_neg = 0;
    curr_true_pos = 0;
    curr_rmse = 0;
    
    start_index = segments(s-1);
    stop_index = segments(s)-1;
    
    curr_ID_labels = ID_labels(start_index:stop_index);
    curr_impact_times = impacts(start_index:stop_index,1)./Fs_fsr;
    
    final_estimates_idx = find(final_estimates_labels(:,1) >= curr_impact_times(1) & final_estimates_labels(:,1) <= curr_impact_times(end));
    curr_estimate_times = final_estimates_labels(final_estimates_idx,1);
    curr_estimate_labels = final_estimates_labels(final_estimates_idx,2);
    
    totalN = length(curr_impact_times);
    for t = 1:totalN
        curr_time = curr_impact_times(t); % real impact time
        curr_label = curr_ID_labels(t);
        diff = abs(curr_estimate_times - curr_time);
        diff_idx = find(diff < 0.05); % 0.05 error threshold!!!!
        same_label_idx = find(curr_estimate_labels == curr_label);
        first_same_label_time = curr_estimate_times(same_label_idx(1));
        last_same_label_time = curr_estimate_times(same_label_idx(end));
        if ~isempty(diff_idx) % if a estimated impact exists close to the real impact time
            est_idx = find(curr_estimate_labels(diff_idx)==curr_label);
            if ~isempty(est_idx) % if the close impact is labeled correctly
                curr_true_pos = curr_true_pos + 1;
                est_times = curr_estimate_times(diff_idx);
                est_time = est_times(est_idx);
                curr_rmse = curr_rmse + (curr_time - est_time)^2;
            else
                % check it's not the beginning and end of the segment
                if first_same_label_time-0.2 < curr_time
                    if last_same_label_time+0.2 > curr_time
                        curr_false_neg = curr_false_neg + 1;
                    end
                end
            end
        else
            if first_same_label_time-0.2 < curr_time
                if last_same_label_time+0.2 > curr_time
                    curr_false_neg = curr_false_neg + 1;
                end
            end
        end
    end
    
    est_totalN = length(curr_estimate_times);
    for t = 1:est_totalN
        curr_est_time = curr_estimate_times(t); % estimated impact time
        curr_est_label = curr_estimate_labels(t);
        diff = abs(curr_impact_times - curr_est_time);
        diff_idx = find(diff < 0.05);
        if isempty(diff_idx)
            same_label_idx = find(curr_ID_labels == curr_est_label);
            first_same_label_time = curr_impact_times(same_label_idx(1));
            last_same_label_time = curr_impact_times(same_label_idx(end));
            % check that it's not the beginning of the segment
            if first_same_label_time-0.2 < curr_est_time
                % check that it's not the beginning of the segment
                if last_same_label_time+0.2 > curr_est_time
                    curr_false_pos = curr_false_pos + 1;
                end
            end
        end
    end
    
    % store
    false_pos(i-1) = curr_false_pos/est_totalN;
    false_neg(i-1) = curr_false_neg/totalN;
    true_pos(i-1) = curr_true_pos/totalN;
    rmse(i-1) = curr_rmse/curr_true_pos;
end



%% perform GMM on these estimated step times
% load processed data 2x workspace files
estimated_scaled_means = zeros(6,4);
scaled_means = zeros(6,4);
real_means = zeros(6,4);
real_std = zeros(6,4);

intervention_arr = ['regular1','limp1','limp2','weight1','weight2','regular2'];
file_root = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Jenny 1\ProcessedData\both_';

for i = 1:6
    load([file_root,intervention_arr(i)])
    for p = 1:2 % for each person
        idx = find(final_estimates_labels(:,2) == i);
        estimated_impact_times = final_estimates_labels(idx,1);
        
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
            estimated_scaled_means(i,person_idx) = mu(1) + (mu(2)-mu(1))*(1- (proportion(1))^2);
            estimated_scaled_means(i,person_idx+1) = mu(2) + (mu(1)-mu(2))*(1- (proportion(2))^2);
        else
            estimated_scaled_means(i,person_idx) = mu(1) + (mu(2)-mu(1))*abs(0.5-proportion(1));
            estimated_scaled_means(i,person_idx+1) = mu(2) + (mu(1)-mu(2))*abs(0.5-proportion(2));
        end
        
        % calculate real step times
        left_right_diff = [];
        right_left_diff = [];
        real_differences = [];
        real_idx = find(impacts(:,4) == person_idx | impacts(:,4) == (person_idx+1));
        real_impact_times = impacts(real_idx,1)./Fs_fsr;
        for i = 2:length(real_impact_times)
            if walk_edges(i) ~= -1 % not the start of a new segment
                real_diff = real_impact_times(i)-real_impact_times(i-1);
                real_differences(end+1) = real_diff;
                if mod(impacts(real_idx(i),4),2) == 0 % right foot
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
            scaled_means(i,person_idx) = mu(1) + (mu(2)-mu(1))*(1- (proportion(1))^2);
            scaled_means(i,person_idx+1) = mu(2) + (mu(1)-mu(2))*(1- (proportion(2))^2);
        else
            scaled_means(i,person_idx) = mu(1) + (mu(2)-mu(1))*abs(0.5-proportion(1));
            scaled_means(i,person_idx+1) = mu(2) + (mu(1)-mu(2))*abs(0.5-proportion(2));
        end
        
        % real step times average & std
        real_means(i,person_idx) = mean(left_right_diff);
        real_means(i,person_idx+1) = mean(right_left_diff);
        real_std(i,person_idx) = std(left_right_diff);
        real_std(i,person_idx+1) = std(right_left_diff);
    end
end

% save(filename,'estimated_scaled_means','scaled_means','real_means','real_std','-append')







































