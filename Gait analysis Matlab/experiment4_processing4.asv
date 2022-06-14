% this file calculates success rate, GMM, create localization csv
%% success rate calculation

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
    
    real_labels = clean_impacts(clean_idx,4);
    real_labels(real_labels == 2) = 1;
    real_labels(real_labels == 3) = 2;
    real_labels(real_labels == 4) = 2;
    real_times = clean_impacts(clean_idx,1)./Fs_fsr;
    
    estimate_idx = find(clean_detected_starts >= start_time & clean_detected_starts <= stop_time);
    estimate_times = clean_detected_starts(estimate_idx);
    estimate_labels = clean_Y_predict(estimate_idx);
    
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

%% perform GMM on these estimated step times for one intervention

estimated_scaled_means = zeros(2,2);
scaled_means = zeros(2,2);
real_means = zeros(2,2);
real_stds = zeros(2,2);
for p = 1:2 % for each person
    idx = find(clean_Y_predict == p);
    estimated_impact_times = clean_detected_starts(idx);

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

%% make localization csv

featureV_p1 = zeros(1,53);
featureV_p2 = zeros(1,53);

% identify the correct estimates
for i = 1:length(clean_impacts(:,1))
    real_time = clean_impacts(i,1)/Fs_fsr;
    real_label = clean_impacts(i,4);
    if real_label == 1 | real_label == 2
        real_label = 1;
    else
        real_label = 2;
    end
    est_idx = find(clean_detected_starts./Fs_pcb > (real_time - 0.05) & clean_detected_starts./Fs_pcb < (real_time + 0.05));
    est_labels = clean_Y_predict(est_idx);
    found_idx = find(est_labels == real_label);
    if ~isempty(found_idx)
        % found a correct estimate
        
        
        
    end
end
    

for i = 1:est_impactN
    curr_time = clean_final_estimates_labels(i,1);
    curr_label = clean_final_estimates_labels(i,2);
    if curr_label == 1 | curr_label == 2
        curr_label = 1;
    else
        curr_label = 2;
    end
    relevant_idx = find(impacts(:,4) == (curr_label-1)*2 + 1 | impacts(:,4) == (curr_label-1)*2 + 2);
    [time_diff,close_minidx] = min(abs(impacts(relevant_idx,1)./Fs_fsr - curr_time));
    if time_diff < 0.05
        correct_estimates(end+1,:) = clean_final_estimates_labels(i,:);
        real_idx = relevant_idx(close_minidx);
        foot_label = impacts(real_idx,4);
        correct_coords(end+1,:) = [coordinates(real_idx,1),coordinates(real_idx,3),foot_label];
        correct_ta(end+1,:) = acc_pks(real_idx,2); % get Y accel (check?)
    end
%     close_arr = find(abs(clean_impacts(relevant_idx,1)./Fs_fsr - curr_time) < 0.05);
%     if ~isempty(close_arr)
%         correct_estimates(end+1,:) = clean_final_estimates_labels(i,:);
%         real_time = clean_impacts(relevant_idx(close_arr(1)),1)/Fs_fsr;
%         foot_label = clean_impacts(relevant_idx(close_arr(1)),4);
%         coords = [allmocap(round(real_time*Fs_mocap),1,foot_label),allmocap(round(real_time*Fs_mocap),3,foot_label),foot_label];
%         correct_coords(end+1,:) = coords;
%     end
end

trainValues = ones(correctN,2); % trainvalue label, person ID
seg_array = unique(correct_estimates(:,3));
% label starts of walking segments
for i = 1:length(segments_list)
    seg = segments_list(i);
    this_seg_idx = find(correct_estimates(:,3) == seg);
    person1_idx = find(correct_estimates(this_seg_idx,2) == 1 | correct_estimates(this_seg_idx,2) == 2);
    person2_idx = find(correct_estimates(this_seg_idx,2) == 3 | correct_estimates(this_seg_idx,2) == 4);
    trainValues(this_seg_idx(person1_idx(1)),:) = [0, 1];
    trainValues(this_seg_idx(person2_idx(1)),:) = [0, 2];
end

% assign k-fold trainValue
segN = seg_array(end)-1;
seg_mod = round(segN/5); % use this as counter for mod calculations
seg_count = 1;
seg_new = false;
for i = 1:correctN
    % if not start of segment
    if trainValues(i,1) ~= 0
        curr_seg = correct_estimates(i,3)-1;
        curr_person = correct_estimates(i,2);
        if curr_person == 1 || curr_person == 2
            curr_person = 1;
        else
            curr_person = 2;
        end
        if mod(curr_seg,seg_mod) ~= 0 | seg_count == 5
            trainValues(i,:) = [seg_count,curr_person];
            seg_new = true;
        else
            if seg_new == true % the first time segment has changed, increment
                seg_count = seg_count + 1;
                seg_new = false;
            end
            trainValues(i,:) = [seg_count,curr_person];
        end
    end
end

for i = 1:correctN
    % take, k-fold trainValue, xcoord, ta, {features}, {feature ratios},
    % prev xcoord
    if trainValues(i,1) == 0
        feature = [1,0,correct_coords(i,1),abs(correct_ta(i)),est_cwt_peak(i,:),est_cwt_energy(i,:),est_peak_mag(i,:),est_energy(i,:),zeros(1,25)];
    else
        % find prev index
        prev_idx = i-1;
        while floor(correct_estimates(prev_idx,2)/3) ~= floor(correct_estimates(i,2)/3)
            prev_idx = prev_idx - 1;
        end
        diff_cwt_peak = est_cwt_peak(i,:) - est_cwt_peak(prev_idx,:);
        diff_cwt_energy = est_cwt_energy(i,:) - est_cwt_energy(prev_idx,:);
        diff_peak_mag = est_peak_mag(i,:) - est_peak_mag(prev_idx,:);
        diff_energy = est_energy(i,:) - est_energy(prev_idx,:);
        feature = [1,trainValues(i,1),correct_coords(i,1),abs(correct_ta(i)),...
            est_cwt_peak(i,:),est_cwt_energy(i,:),est_peak_mag(i,:),est_energy(i,:),...
            diff_cwt_peak,diff_cwt_energy,diff_peak_mag,diff_energy,correct_coords(prev_idx,1)];
    end
    if correct_estimates(i,2) == 1 | correct_estimates(i,2) == 2
        featureV_p1(end+1,:) = feature;
    else
        featureV_p2(end+1,:) = feature;
    end
end
featureV_p1(1,:) = [];
featureV_p2(1,:) = [];
filename1 = [filename,'_localization_p1_withta.csv'];
filename2 = [filename,'_localization_p2_withta.csv'];
% filename1 = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Jenny 1\ProcessedData\ExcelData\both_regular1_localization_p1_withta.csv';
% filename2 = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Jenny 1\ProcessedData\ExcelData\both_regular1_localization_p2_withta.csv';
writematrix(featureV_p1,filename1)
writematrix(featureV_p2,filename2)
