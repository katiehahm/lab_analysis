% this file calculates success rate, GMM, create localization csv
%% success rate calculation IGNORE?

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

processedfilepath = 'C:\Users\katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Jenny 1\ProcessedData\both_limp2.mat';
load(processedfilepath)

scaled_means = zeros(2,2);
real_means = zeros(2,2);
real_stds = zeros(2,2);
for p = 1:2 % for each person
    idx = find(clean_est_impacts(:,2) == p);
    estimated_impact_times = clean_est_impacts(idx,1)./Fs_pcb;

    % calculate all step times from impact times
    step_times = [];
    for x = 2:length(estimated_impact_times)
        curr_diff = estimated_impact_times(x) - estimated_impact_times(x-1);
        if curr_diff < 1.5 % not a turning point
            step_times(end+1) = curr_diff;
        end
    end
    
    if p == 1
        filename = [processedfilepath(1:end-4),'_p1_GMM.csv'];
        writematrix(step_times.',filename)
    else
        filename = [processedfilepath(1:end-4),'_p2_GMM.csv'];
        writematrix(step_times.',filename)
    end

    % GM analysis from estimated step times
    GM = fitgmdist(transpose(step_times),2,'RegularizationValue',0.000001);
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
    estimated_scaled_means(p,:) = [min(mean1,mean2), max(mean1,mean2)];
        
    % calculate real step times
    left_right_diff = [];
    right_left_diff = [];
    real_differences = [];
    real_idx = find(floor(clean_impacts(:,2)/10) == p);
    real_impact_times = clean_impacts(real_idx,1);
    for x = 2:length(real_impact_times)
        real_diff = real_impact_times(x)-real_impact_times(x-1);
        if real_diff < 1.5 % not a turning point
            real_differences(end+1) = real_diff;
            if mod(clean_impacts(real_idx(x),2),10) == 2 % right foot
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

scaled_means
real_means

% then run python GMM_derivation.py

%% from python GMM
estimated_scaled_means = [0.5262,0.7829;0.5460,0.5983];
save(processedfilepath,'processedfilepath','estimated_scaled_means','scaled_means','real_means','real_stds','-append')

%% make localization csv

processedfilepath = 'C:\Users\katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Jenny 1\ProcessedData\both_limp2.mat';
load(processedfilepath)
featureV_p1 = zeros(1,32);
featureV_p2 = zeros(1,32);
time_buffer = 0.05;

% identify the correct estimates
for i = 1:length(clean_impacts(:,1))
    real_time = clean_impacts(i,1);
    real_label = clean_impacts(i,2);
    real_label = floor(real_label/10);
    est_idx = find(clean_est_impacts(:,1)./Fs_pcb > (real_time - time_buffer) & clean_est_impacts(:,1)./Fs_pcb < (real_time + time_buffer));
    est_labels = clean_est_impacts(est_idx,2);
    found_idx = find(est_labels == real_label);
    if ~isempty(found_idx)
        % found a correct estimate
        new_array = [clean_impacts(i,3),clean_impacts(i,5),clean_est_impacts(est_idx(found_idx),15:44)];
        if real_label == 1
            featureV_p1 = [featureV_p1; new_array];
        else
            featureV_p2 = [featureV_p2; new_array];
        end
    end
end
featureV_p1(1,:) = [];
featureV_p2(1,:) = [];

filename1 = [processedfilepath(1:end-4),'_localization_p1_withta.csv'];
filename2 = [processedfilepath(1:end-4),'_localization_p2_withta.csv'];
writematrix(featureV_p1,filename1)
writematrix(featureV_p2,filename2)
