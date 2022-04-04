%% used to detect footfalls using cwt
% 3/24/22

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Praneeth experiment 3_11_22\ProcessedData\both_towards1';
load(string(data_root_katie))

freq_lower = 100; % Hz; frequency limits for cwt
freq_higher = 350;

final_estimates = [0,0]; % [impact time, segment ID num]
impact_thresh = 0.04; % (s) of how close an impact should be to the last

arrival_window = 3000; % 5000/Fs_pcb width of window to find arrival time
impact_time_thresh = 0.1; % if two impacts are 0.1s within each other, combine
correct_thresh = 0.05; % if predicted is within this thresh to real impact time, it's correct

for w = 2:length(segments) % for each walking segment
    start_time = impacts(segments(w-1)+1,1)/Fs_fsr-1;
    stop_time = impacts(segments(w),1)/Fs_fsr+0.1;
    
    start_idx_pcb = findTindex(start_time,pcbTime);
    stop_idx_pcb = findTindex(stop_time,pcbTime);
    
    impacttimes = impacts(segments(w-1)+1:segments(w),1)/Fs_fsr;
    
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
        impacttimes = impacts(segments(w-1)+1:segments(w),1)/Fs_fsr;
        plot((impacttimes-start_time)*Fs_pcb,0,'rx','MarkerSize',8)
        for i = 1:length(locs)
            if i == 1
                starti = max(locs(i) - arrival_window, 1);
            else
                [minval, minidx] = min(sum_smooth_cwt(locs(i-1):locs(i)));
                starti = max(locs(i) - arrival_window, locs(i-1)+minidx - arrival_window/5);
            end
            if minval > 0.001 && i ~= 1
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

%% using above, look at all combinations, find one with min uncertainty 3/30/22

close all
for s = 3:3%length(segments)
    start_index = segments(s-1)+1;
    stop_index = segments(s);

    step_times_idx = find(final_estimates(:,2) == s);
    step_times = final_estimates(step_times_idx,1)./Fs_pcb;
    
    curr_ID_labels = ID_labels(start_index:stop_index);
    real_impact_times = impacts(start_index:stop_index,1)./Fs_fsr;
    
    [estimateID, step_times] = recursive_stepID(curr_ID_labels, real_impact_times, sort(step_times),step_o,step_x);
    
    title(sprintf('Final estimated ID labels %d segment'),s)
    
    
    % calculating success rate
    successN = length(real_impact_times);
    success_count = 0;
    rmse_val = 0;
    for i = 1:successN
        real_ID = curr_ID_labels(i);
        real_time = real_impact_times(i);
        est_thispp = find(estimateID == real_ID);
        diff_array = abs(step_times(est_thispp)-real_time);
        if ~isempty(find(diff_array < 0.15))
            success_count = success_count + 1;
            [minval, ~] = min(diff_array);
            rmse_val = rmse_val + minval;
        end
    end
    rmse = rmse_val/successN
    success = success_count/successN
    extra_detect = length(step_times) - length(real_impact_times)
    
end

seg3_estimated_steptimes = step_times;
seg3_estimateID = estimateID;

%% save estimated step ID date 4/4/22

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Praneeth experiment 3_11_22\ProcessedData\both_towards1_estimatedStepTimes';
save(data_root_katie,'seg2_estimateID', 'seg2_estimated_steptimes',...
    'seg3_estimateID', 'seg3_estimated_steptimes',...
    'seg4_estimateID', 'seg4_estimated_steptimes',...
    'seg5_estimateID', 'seg5_estimated_steptimes',...
    'seg6_estimateID', 'seg6_estimated_steptimes',...
    'seg7_estimateID', 'seg7_estimated_steptimes',...
    'seg8_estimateID', 'seg8_estimated_steptimes',...
    'seg9_estimateID', 'seg9_estimated_steptimes',...
    'seg10_estimateID', 'seg10_estimated_steptimes',...
    'seg11_estimateID', 'seg11_estimated_steptimes')

%% perform GMM on these estimated step times




