% for use after experiment4_allprocessingcompiled.m
% gets features and makes csv for localization estimation
%% feature extraction 5/11/22

% constants
est_impactN = length(clean_final_estimates_labels);
window_thresh = 0.0825; % seconds of how long an impact lasts
start_idx_thresh = 1000; % count 1000 idx before curr_time to do aicpick
freq_lower = 100; % Hz; frequency limits for cwt
freq_higher = 350;
noise_thresh = [0.000285,0.00026,0.00018,0.00026,0.00024,0.00056]; % determined empirically by looking at plot

correct_estimates = [0,0,0];
% find all correct impacts
for i = 1:est_impactN
    curr_time = clean_final_estimates_labels(i,1);
    curr_label = clean_final_estimates_labels(i,2);
    relevant_idx = find(clean_impacts(:,4) == (curr_label-1)*2 + 1 | clean_impacts(:,4) == (curr_label-1)*2 + 2);
    if ~isempty(find(abs(clean_impacts(relevant_idx,1)./Fs_fsr - curr_time) < 0.05))
        correct_estimates(end+1,:) = clean_final_estimates_labels(i,:);
    end
end
correct_estimates(1,:) = [];
correctN = length(correct_estimates);

% initialize feature arrays
est_coordinates = zeros(correctN,2);
est_arrival_idx = zeros(correctN,sensorN);
est_last_idx = zeros(correctN,sensorN);
est_peak_idx = zeros(correctN,sensorN);
est_peak_mag = zeros(correctN,sensorN);
est_energy = zeros(correctN,sensorN);
est_cwt_energy = zeros(correctN,sensorN);
est_overlapping = zeros(correctN,1); % value = 1 if two impacts overlap

% filter pcb signal
filt_pcbD = lpf_data(pcbData);

% extract values
for i = 1:correctN
    curr_time = correct_estimates(i,1);
    curr_idx = round(curr_time*Fs_pcb);
    % check if overlapping
    if i~= correctN & curr_time == correct_estimates(i+1,1)
        est_overlapping(i) = 1;
        est_overlapping(i+1) = 1;
    % if not overlapping
    elseif est_overlapping(i) == 0
        for s = 1:sensorN
            % extract features from this impact
            curr_label = correct_estimates(i,2);
            % define window end
            if i ~= correctN % if it isn't the last impact
                window_end = correct_estimates(i+1,1);
            else % last impact
                window_end = length(pcbData)/Fs_pcb;
            end
            window = filt_pcbD(curr_idx:round(window_end*Fs_pcb),s);
            noise = noise_thresh(s);
            if ~isempty(window)
                [up,lo] = envelope(window,300,'peak');
                indeces = find(up < noise);
                while isempty(indeces)
                    noise = noise + 0.00001;
                    indeces = find(up < noise);
                end
                last_idx = indeces(1) + curr_idx;
            end
            est_last_idx(i,s) = last_idx;
            % define window start
            start_idx = round(curr_idx - start_idx_thresh);
            curr_window = filt_pcbD(start_idx:last_idx,s);
            % extract arrival time
            est_arrival_idx(i,s) = aic_pick(curr_window, 'to_peak')+start_idx;
            curr_window = filt_pcbD(est_arrival_idx(i,s):last_idx,s);
            % extract peak mag
            [maxval, maxidx] = max(curr_window);
            est_peak_mag(i,s) = maxval;
            est_peak_idx(i,s) = maxidx + start_idx - 1;
            % extract energy
            est_energy(i,s) = sum(abs(curr_window).^2);
            % extract cwt energy
            [wt,f] = cwt(curr_window, Fs_pcb);
            valid_f_idx = find(freq_higher & f > freq_lower);
            cwt_mag = abs(wt(valid_f_idx,:));
            sum_cwt = sum(cwt_mag,1);
            sum_smooth_cwt = movmean(sum_cwt, 800);
            est_cwt_energy(i,s) = sum(abs(sum_smooth_cwt));
        end
        % extract ground truth coordinates
        Limpact_idx = find(clean_impacts(:,4) == ((curr_label-1)*2 + 1));
        Rimpact_idx = find(clean_impacts(:,4) == ((curr_label-1)*2 + 2));
        [Lminval,Lcorresponding_real_idx] = min(abs(clean_impacts(Limpact_idx,1)/Fs_fsr - curr_time));
        [Rminval,Rcorresponding_real_idx] = min(abs(clean_impacts(Rimpact_idx,1)/Fs_fsr - curr_time));
        mocap_idx = round(Fs_mocap*(curr_time + 0.05)); % delay for the heel to get "settled"
        if Lminval < Rminval
            mocap_label = (curr_label - 1)*2 + 2;
        else
            mocap_label = (curr_label - 1)*2 + 1;
        end
        est_coordinates(i,:) = [allmocap(mocap_idx,1,mocap_label), allmocap(mocap_idx,3,mocap_label)];
    end
end

%% trying above feature extraction w/ cwt approach, similar code 5/14/22

% constants
est_impactN = length(clean_final_estimates_labels);
window_thresh = 0.0825; % seconds of how long an impact lasts
start_idx_thresh = 1000; % count 1000 idx before curr_time to do aicpick
freq_lower = 100; % Hz; frequency limits for cwt
freq_higher = 350;
noise_thresh = [0.000285,0.00026,0.00018,0.00026,0.00024,0.00056]; % determined empirically by looking at plot

correct_estimates = [0,0,0];
correct_coords = [0,0,0];
% find all correct impacts & ground truth coords
for i = 1:est_impactN
    curr_time = clean_final_estimates_labels(i,1);
    curr_label = clean_final_estimates_labels(i,2);
    relevant_idx = find(clean_impacts(:,4) == (curr_label-1)*2 + 1 | clean_impacts(:,4) == (curr_label-1)*2 + 2);
    close_arr = find(abs(clean_impacts(relevant_idx,1)./Fs_fsr - curr_time) < 0.05);
    if ~isempty(close_arr)
        correct_estimates(end+1,:) = clean_final_estimates_labels(i,:);
        real_time = clean_impacts(relevant_idx(close_arr(1)),1)/Fs_fsr;
        foot_label = clean_impacts(relevant_idx(close_arr(1)),4);
        coords = [allmocap(round(real_time*Fs_mocap),1,foot_label),allmocap(round(real_time*Fs_mocap),3,foot_label),foot_label];
        correct_coords(end+1,:) = coords;
    end
end
correct_estimates(1,:) = [];
correct_coords(1,:) = [];
correctN = length(correct_estimates);

% initialize feature arrays
est_coordinates = zeros(correctN,3);
est_arrival_idx = zeros(correctN,sensorN);
est_last_idx = zeros(correctN,sensorN);
est_peak_idx = zeros(correctN,sensorN);
est_peak_mag = zeros(correctN,sensorN);
est_energy = zeros(correctN,sensorN);
est_cwt_energy = zeros(correctN,sensorN);
est_cwt_peak = zeros(correctN,sensorN);
est_overlapping = zeros(correctN,1); % value = 1 if two impacts overlap

% filter pcb signal
filt_pcbD = lpf_data(pcbData);

segments_list = unique(correct_estimates(:,3));
segmentsN = segments_list(end)-segments_list(1) + 1;
% extract values
for seg = 1:segmentsN
    seg_idx = find(correct_estimates(:,3) == (seg+1));
    seg_estimates = correct_estimates(seg_idx,:);
    estimatesN = length(seg_estimates(:,1));
    for s = 1:sensorN
        % extract cwt of segment
        not_overlap = [];
        first_idx = round((seg_estimates(1,1)-0.5)*Fs_pcb);
        last_idx = round((seg_estimates(estimatesN,1)+0.5)*Fs_pcb);
        seg_pcb = filt_pcbD(first_idx:last_idx,s);
        [wt,f] = cwt(seg_pcb,Fs_pcb);
        valid_f_idx = find(freq_higher & f > freq_lower);
        cwt_mag = abs(wt(valid_f_idx,:));
        sum_cwt = sum(cwt_mag,1);
        sum_smooth_cwt = movmean(sum_cwt, 800);
        pcb2cwt_idx = first_idx; % add this to result to get pcb idx
        % extract features from cwt sum
        for i = 1:estimatesN
            % check for overlapping
            if i~= estimatesN & seg_estimates(i,1) == seg_estimates(i+1,1)
                est_overlapping(seg_idx(i)) = 1;
                est_overlapping(seg_idx(i+1)) = 1;
            elseif est_overlapping(seg_idx(i)) == 0
                not_overlap(end+1) = i;
                arrival_idx = round(seg_estimates(i,1)*Fs_pcb);
%                 starti = round((seg_estimates(i,1)-0.035)*Fs_pcb);
                if i ~= estimatesN
                    [lasti,nextimpact] = min([round(seg_estimates(i+1,1)*Fs_pcb),arrival_idx + round(0.46*Fs_pcb)]);
                else
                    lasti = min([last_idx,arrival_idx + round(0.31*Fs_pcb)]);
                    nextimpact = 2;
                end
                window = sum_smooth_cwt(arrival_idx-pcb2cwt_idx+1:lasti-pcb2cwt_idx+1);
                deriv_window = movmean(diff(window),1200);
                [~,minidx] = min(deriv_window);
                [~,maxidx] = max(deriv_window);
                % arrival idx is the same for all impacts
%                 arrival_idx = round(seg_estimates(i,1)*Fs_pcb);
                est_arrival_idx(seg_idx(i),s) = arrival_idx;
                % window end time
                [~,zeroidx] = min(abs(deriv_window(minidx:end)));
                abs_last_idx = zeroidx + minidx - 1 + arrival_idx;

%                 zerothresh = 0.000001;
%                 zeroidx = find(abs(deriv_window(minidx:end)) < zerothresh);
%                 while isempty(zeroidx)
%                     zerothresh = zerothresh + 0.000001;
%                     zeroidx = find(abs(deriv_window(minidx:end)) < zerothresh);
%                 end
%                 abs_last_idx = zeroidx(1) + minidx - 1 + arrival_idx;

                % flip the window and do aic pick to find last impact
    %             flip_window = flip(window);
    %             flip_last_idx = aic_pick(flip_window, 'to_peak');
    %             abs_last_idx = length(window) - flip_last_idx + starti;
                est_last_idx(seg_idx(i),s) = abs_last_idx;
                % cwt energy
                if (abs_last_idx - arrival_idx) > 10 % sometimes peak is so small it's a straight line, then last idx is right after arrival idx
                    window2 = filt_pcbD(arrival_idx:abs_last_idx,s);
                    [wt2,f2] = cwt(window2,Fs_pcb);
                    valid_f_idx2 = find(freq_higher & f2 > freq_lower);
                    cwt_mag2 = abs(wt2(valid_f_idx2,:));
                    sum_cwt2 = sum(cwt_mag2,1);
                    sum_smooth_cwt2 = movmean(sum_cwt2, 800);
                    est_cwt_energy(seg_idx(i),s) = sum(sum_smooth_cwt2);
                    % peak mag
                    est_peak_mag(seg_idx(i),s) = max(window);
                    % peak cwt mag
                    est_cwt_peak(seg_idx(i),s) = max(sum_smooth_cwt2);
                    % energy
                    est_energy(seg_idx(i),s) = sum(abs(window).^2);
                end
                % ground truth coordinates
                curr_time = correct_estimates(seg_idx(i),1);
                curr_label = correct_estimates(seg_idx(i),2);
                Limpact_idx = find(clean_impacts(:,4) == ((curr_label-1)*2 + 1));
                Rimpact_idx = find(clean_impacts(:,4) == ((curr_label-1)*2 + 2));
                [Lminval,Lcorresponding_real_idx] = min(abs(clean_impacts(Limpact_idx,1)/Fs_fsr - curr_time));
                [Rminval,Rcorresponding_real_idx] = min(abs(clean_impacts(Rimpact_idx,1)/Fs_fsr - curr_time));
                mocap_idx = round(Fs_mocap*(curr_time + 0.05)); % delay for the heel to get "settled"
                if Lminval < Rminval
                    mocap_label = (curr_label - 1)*2 + 2;
                else
                    mocap_label = (curr_label - 1)*2 + 1;
                end
                est_coordinates(i,:) = [allmocap(mocap_idx,1,mocap_label), allmocap(mocap_idx,3,mocap_label),mocap_label];
            end
        end
        % sanity check per sensor
%         figure;
%         plot(sum_smooth_cwt)
%         hold on
%         plot(correct_estimates(seg_idx(not_overlap),1)*Fs_pcb - first_idx,sum_smooth_cwt(round(correct_estimates(seg_idx(not_overlap),1)*Fs_pcb - first_idx)),'gx','MarkerSize',14)
%         plot(est_arrival_idx(seg_idx(not_overlap),s)-first_idx,sum_smooth_cwt(round(est_arrival_idx(seg_idx(not_overlap),s)-first_idx)),'rx','MarkerSize',11)
%         plot(est_last_idx(seg_idx(not_overlap),s)-first_idx,sum_smooth_cwt(round(est_last_idx(seg_idx(not_overlap),s)-first_idx)),'kx','MarkerSize',8)
%         figure; plot(movmean(diff(sum_smooth_cwt),1200))
    end
end

%% sanity check 5/13/22

% to check if overlapping detection is correct
% red x's should be vertically stacked
figure;
plot(correct_estimates(:,1),linspace(1,correctN,correctN),'b.')
hold on
plot(correct_estimates(find(est_overlapping == 1),1),find(est_overlapping == 1),'rx')

% to test that missing features are only due to overlapping impacts
% should run bc the arrays are same length
arrival_zeros = find(est_arrival_idx(:,1) == 0);
figure;
plot(arrival_zeros,find(est_overlapping == 1))

% to check if impact window is accurate
arrival_nonzeros = find(est_arrival_idx(:,1) ~= 0);
for s = 1:sensorN
    figure;
    plot(pcbTime, filt_pcbD(:,s))
    hold on
    plot(est_arrival_idx(arrival_nonzeros,s)./Fs_pcb,0,'rx','MarkerSize',14)
    plot(est_last_idx(arrival_nonzeros,s)./Fs_pcb,0,'gx','MarkerSize',9)
%     plot(correct_estimates(:,1),0,'cx')
end

% check if coordinates are correct
for i = 1:4
    figure;
    plot(allmocap(:,1,i)) % foot i
    hold on
    foot_idx = find(correct_coords(:,3) == i);
    est_coord = correct_coords(foot_idx,1);
    est_arrive = est_arrival_idx(foot_idx,1)/Fs_pcb;
    plot(est_arrive*Fs_mocap,est_coord,'rx')
end


%% wiener filter & interpolation for overlapping impacts 5/12/22

overlap_idx = find(est_overlapping == 1);
for i = 7:length(overlap_idx) % put 7 here bc first segment is trash
    curr_idx = overlap_idx(i); % idx within length of clean_final...
    data_idx = curr_idx*Fs_pcb; % idx within data segment
    curr_time = clean_final_estimates_labels(curr_idx,1);
    curr_label = clean_final_estimates_labels(curr_idx,2);
    prev_labels = clean_final_estimates_labels(curr_idx-4:curr_idx-1,2);
    prev_idx = find(prev_labels == curr_label);
    % last_idx is idx of impacts, not idx of data
    last_idx = curr_idx - 4 + prev_idx(end) - 1;
    for s = 1:sensorN
        % previous clip
        prev_clip = filt_pcbD(est_arrival_idx(last_idx,s):est_last_idx(last_idx,s),s);
        % get overlapping clip
        window = filt_pcbD(data_idx:data_idx + Fs_pcb/2,s); % set window of 0.5s
        noise = noise_thresh(s);
        if ~isempty(window)
            [up,lo] = envelope(window,300,'peak');
            indeces = find(up < noise);
            while isempty(indeces)
                noise = noise + 0.00001;
                indeces = find(up < noise);
            end
            data_last_idx = indeces(1) + data_idx;
        end
        start_idx = round(data_idx - start_idx_thresh);
        curr_window = filt_pcbD(start_idx:data_last_idx,s);
        arrival_idx = aic_pick(curr_window, 'to_peak')+start_idx;
        window = filt_pcbD(arrival_idx:data_last_idx,s);
        est_arrival_idx(curr_idx,s) = arrival_idx;
        est_last_idx(curr_idx,s) = data_last_idx;
        % make clips same length
        prev_clipL = length(prev_clip);
        clipL = length(window);
        if prev_clipL > clipL % if previous clip is longer than overlapping clip
            % extend current clip
            L_diff = prev_clipL - clipL;
            window = filt_pcbD(arrival_idx - round(L_diff/2):data_last_idx + (L_diff - round(L_diff/2)),s);
        else % previous clip is shorter than current overlapping clip
            L_diff = clipL - prev_clipL;
            prev_clip = filt_pcbD(est_arrival_idx(last_idx,s) - round(L_diff/2):est_last_idx(last_idx,s)+ (L_diff - round(L_diff/2)),s);
            % set any added parts above noise threshold to be noise
            make_noise = find(prev_clip(1:round(L_diff/2)) > noise_thresh(s));
            prev_clip(make_noise) = noise_thresh(s);
            make_noise = find(prev_clip(est_last_idx(last_idx,s):end) > noise_thresh(s));
            prev_clip(make_noise + est_last_idx(last_idx,s) - 1) = noise_thresh(s);
        end
        % perform wiener filter
        [yhat, ~] = wienerFilter(prev_clip,window,1,true,Fs_pcb);
        % interpolate to find scaling factor
        next_labels = clean_final_estimates_labels(curr_idx+1:curr_idx+4,2);
        next_idx_list = find(next_labels == curr_label);
        next_idx = curr_idx + 1 + next_idx_list(1) - 1;
        next_peak_mag = est_peak_mag(next_idx,s);
        prev_peak_mag = est_peak_mag(last_idx,s);
        curr_peak_mag = (next_peak_mag+prev_peak_mag)/2;
        est_peak_mag(curr_idx,s) = curr_peak_mag;
        scaling_factor = curr_peak_mag/max(yhat);
        yhat = yhat.*scaling_factor;
        % extract energy from extracted signal
        est_energy(curr_idx,s) = sum(abs(yhat).^2);
        % extract cwt energy from extracted signal
        [wt,f] = cwt(yhat, Fs_pcb);
        valid_f_idx = find(freq_higher & f > freq_lower);
        cwt_mag = abs(wt(valid_f_idx,:));
        sum_cwt = sum(cwt_mag,1);
        sum_smooth_cwt = movmean(sum_cwt, 800);
        est_cwt_energy(curr_idx,s) = sum(abs(sum_smooth_cwt));
    end
end

































