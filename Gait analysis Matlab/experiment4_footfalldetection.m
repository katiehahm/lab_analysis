%% used to detect footfalls using cwt
% 3/24/22

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Praneeth experiment 3_11_22\ProcessedData\both_towards1';
load(string(data_root_katie))

% find walking segments


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
% load processed data 2x workspace files
GMmodels = zeros(2,9);
scaled_means = zeros(2,3);

for i = 1:2 % for each person
    differences = [];
    
    for s = 2:11
        if s == 2
            idx = find(seg2_estimateID == i);
            steptimes = seg2_estimated_steptimes(idx);
        elseif s == 3
            idx = find(seg3_estimateID == i);
            steptimes = seg3_estimated_steptimes(idx);
        elseif s == 4
            idx = find(seg4_estimateID == i);
            steptimes = seg4_estimated_steptimes(idx);
        elseif s == 5
            idx = find(seg5_estimateID == i);
            steptimes = seg5_estimated_steptimes(idx);
        elseif s == 6
            idx = find(seg6_estimateID == i);
            steptimes = seg6_estimated_steptimes(idx);
        elseif s == 7
            idx = find(seg7_estimateID == i);
            steptimes = seg7_estimated_steptimes(idx);
        elseif s == 8
            idx = find(seg8_estimateID == i);
            steptimes = seg8_estimated_steptimes(idx);
        elseif s == 9
            idx = find(seg9_estimateID == i);
            steptimes = seg9_estimated_steptimes(idx);
        elseif s == 10
            idx = find(seg10_estimateID == i);
            steptimes = seg10_estimated_steptimes(idx);
        elseif s == 11
            idx = find(seg11_estimateID == i);
            steptimes = seg11_estimated_steptimes(idx);
        end
            
        for t = 2:length(steptimes)
            differences(end+1) = steptimes(t) - steptimes(t-1);
        end
    end

    % estimated step times GMM results
    GM = fitgmdist(transpose(differences),2,'RegularizationValue',0.000001);
    proportion = GM.ComponentProportion;
    mu = GM.mu;
    sig = GM.Sigma;
    GMmodels(i,:) = [mu(1), mu(2), sig(1), sig(2), proportion(1), proportion(2), GM.AIC, GM.BIC, GM.NegativeLogLikelihood];
    
    [prop1,~] = max(GMmodels(i,5:6));
    if (1-prop1) < abs(0.5-prop1)
        scaled_means(i,1) = GMmodels(i,1) + (GMmodels(i,2)-GMmodels(i,1))*(1- (GMmodels(i,5))^2);
        scaled_means(i,2) = GMmodels(i,2) + (GMmodels(i,1)-GMmodels(i,2))*(1- (GMmodels(i,6))^2);
        scaled_means(i,3) = abs(scaled_means(i,1)-scaled_means(i,2));
    else
        scaled_means(i,1) = GMmodels(i,1) + (GMmodels(i,2)-GMmodels(i,1))*abs(0.5- GMmodels(i,5));
        scaled_means(i,2) = GMmodels(i,2) + (GMmodels(i,1)-GMmodels(i,2))*abs(0.5- GMmodels(i,6));
        scaled_means(i,3) = abs(scaled_means(i,1)-scaled_means(i,2));
    end
    i
    disp("Estimated scaled means")
    scaled_means
    
    % real left and right average step times
    left_right_diff = [];
    right_left_diff = [];
    real_differences = [];
    for s = 2:length(segments)
        start_index = segments(s-1)+1;
        stop_index = segments(s);

        real_impact_time_arr = impacts(start_index:stop_index,1)./Fs_fsr;
        [real_impact_time_arr,sort_idx] = sort(real_impact_time_arr);
        impact_ID = impacts(start_index:stop_index,4);
        impact_ID = impact_ID(sort_idx);
        this_impact_ID_idx = find(floor(impact_ID./3) == (i-1));
        
        for t = 2:length(this_impact_ID_idx)
            curr_i = this_impact_ID_idx(t);
            past_i = this_impact_ID_idx(t-1);
            real_differences(end+1) = (real_impact_time_arr(curr_i) - real_impact_time_arr(past_i));
            if mod(impact_ID(curr_i),2) == 1 % left foot
                left_right_diff(end+1) = (real_impact_time_arr(curr_i) - real_impact_time_arr(past_i));
            else % right foot
                right_left_diff(end+1) = (real_impact_time_arr(curr_i) - real_impact_time_arr(past_i));
            end
        end
        
        
    end
    left_to_right_groundtruth = mean(left_right_diff)
    right_to_left_groundtruth = mean(right_left_diff)
    
    % GMM results of real step times
    mu = mean(real_differences);
    sig = std(real_differences);
    mu_s = [mu; mu+0.01];
    sigma_s = zeros(1,1,2);
    sigma_s(1,1,:) = [sig; sig];
    pcomponents = [1/2,1/2];
    
    GM = fitgmdist(transpose(real_differences),2,'RegularizationValue',0.000001);
    proportion = GM.ComponentProportion;
    mu = GM.mu;
    sig = GM.Sigma;
    GMmodel_real = [mu(1), mu(2), sig(1), sig(2), proportion(1), proportion(2), GM.AIC, GM.BIC, GM.NegativeLogLikelihood];
    
    [prop1,~] = max(GMmodel_real(5:6));
    if (1-prop1) < abs(0.5-prop1)
        scaled_mean1 = GMmodel_real(1) + (GMmodel_real(2)-GMmodel_real(1))*(1- (GMmodel_real(5))^2);
        scaled_mean2 = GMmodel_real(2) + (GMmodel_real(1)-GMmodel_real(2))*(1- (GMmodel_real(6))^2);
        scaled_mean_diff = abs(scaled_mean1-scaled_mean2);
    else
        scaled_mean1 = GMmodel_real(1) + (GMmodel_real(2)-GMmodel_real(1))*abs(0.5- GMmodel_real(5));
        scaled_mean2 = GMmodel_real(2) + (GMmodel_real(1)-GMmodel_real(2))*abs(0.5- GMmodel_real(6));
        scaled_mean_diff = abs(scaled_mean1-scaled_mean2);
    end
    disp("Real step times GMM")
    [scaled_mean1, scaled_mean2, scaled_mean_diff]
end

%% save estimated_impacts (4/5/22)

estimated_impacts = zeros(1, 2);

for s = 2:11
    if s == 2
        estimateID = seg2_estimateID;
        steptimes = seg2_estimated_steptimes;
    elseif s == 3
        estimateID = seg3_estimateID;
        steptimes = seg3_estimated_steptimes;
    elseif s == 4
        estimateID = seg4_estimateID;
        steptimes = seg4_estimated_steptimes;
    elseif s == 5
        estimateID = seg5_estimateID;
        steptimes = seg5_estimated_steptimes;
    elseif s == 6
        estimateID = seg6_estimateID;
        steptimes = seg6_estimated_steptimes;
    elseif s == 7
        estimateID = seg7_estimateID;
        steptimes = seg7_estimated_steptimes;
    elseif s == 8
        estimateID = seg8_estimateID;
        steptimes = seg8_estimated_steptimes;
    elseif s == 9
        estimateID = seg9_estimateID;
        steptimes = seg9_estimated_steptimes;
    elseif s == 10
        estimateID = seg10_estimateID;
        steptimes = seg10_estimated_steptimes;
    elseif s == 11
        estimateID = seg11_estimateID;
        steptimes = seg11_estimated_steptimes;
    end

    for t = 1:length(estimateID)
        add_array = [steptimes(t),estimateID(t)];
        estimated_impacts = [estimated_impacts; add_array];
    end
end
estimated_impacts(1,:) = []; % from initialization

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Praneeth experiment 3_11_22\ProcessedData\both_towards1_estimatedStepTimes';
save(data_root_katie,'estimated_impacts','-append')










































