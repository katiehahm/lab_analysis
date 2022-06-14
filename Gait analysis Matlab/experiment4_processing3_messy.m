%% use decision tree to classify impacts
% load training data
load('C:\Users\katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\both_regular1.mat')
X_train = X;
Y_train = Y;
% load testing data
load('C:\Users\katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\both_regular2.mat')

% detect impact starts
pkprom = [0.0015,0.0015,0.0015,0.0015,0.0015,0.0015];
pkheight = [0.006,0.016,0.0014,0.004,0.0127,0.009];
freq_lower = 100; % Hz; frequency limits for cwt
freq_higher = 350;
arrival_window = 1500; % 5000/Fs_pcb width of window to find arrival time
impact_time_thresh = 0.05; % if two impacts are 0.05s within each other, combine
correct_thresh = 0.05; % if predicted is within this thresh to real impact time, it's correct
estimated_impacts = [0,0,0]; % [impact time, sensor N, segmentN

for seg = 2:length(segments)
    start_time = impacts(segments(seg-1),1)-0.15;
    stop_time = impacts(segments(seg)-1,1)+0.15;

    start_idx_pcb = findTindex(start_time,pcbTime);
    stop_idx_pcb = findTindex(stop_time,pcbTime);

    impacttimes = impacts(segments(seg-1):segments(seg)-1,1);
    
    for s = 1:6 % for each sensor
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
        impacttimes = impacts(segments(seg-1):segments(seg)-1,1);
        o_impactidx = find(floor(impacts(segments(seg-1):segments(seg)-1,2)/10) == 1);
        x_impactidx = find(floor(impacts(segments(seg-1):segments(seg)-1,2)/10) == 2);
%         o_impactidx = find(impacts(segments(seg-1):segments(seg)-1,4) == 1 | impacts(segments(seg-1):segments(seg)-1,4) == 2);
%         x_impactidx = find(impacts(segments(seg-1):segments(seg)-1,4) == 3 | impacts(segments(seg-1):segments(seg)-1,4) == 4);
        plot((impacttimes(o_impactidx)-start_time)*Fs_pcb,0,'ro','MarkerSize',8)
        plot((impacttimes(x_impactidx)-start_time)*Fs_pcb,0,'rx','MarkerSize',8)
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
            % if window start is higher than detected peak, move start up
            while sum_smooth_cwt(starti) > sum_smooth_cwt(locs(i))
                starti = starti + 50;
            end
            window = sum_smooth_cwt(starti:locs(i));
            arrive_idx = aic_pick(window, 'to_peak')+starti - 1;
            add_array = [arrive_idx + start_idx_pcb - 1, s,seg];
            estimated_impacts = [estimated_impacts; add_array];
            plot(arrive_idx,0,'bx') 
        end
    end
end
estimated_impacts(1,:) = []; % from initialization

%% sanity check detected impacts
figure;
plot(impacts(:,1),0,'rx','MarkerSize',8)
hold on
plot(estimated_impacts(:,1)./Fs_pcb,0.5,'bx')
ylim([-0.5 1])
legend('real impacts','estimated impacts')

%% combine impacts
estimates = sort(estimated_impacts(:,1));
final_pred = [];
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
    final_pred(end+1) = mean(estimates(groupidx));
    scores(end+1) = length(groupidx);
end

delete_pred = [];
for i = 2:length(final_pred)
    if final_pred(i) - final_pred(i-1) < thresh
        final_pred(i) = mean(final_pred(i-1:i));
        delete_pred(end+1) = i-1;
        scores(i) = scores(i) + scores(i-1);
    end
end
final_pred(delete_pred) = [];
scores(delete_pred) = [];

figure;
plot(final_pred./Fs_pcb,0.5,'bx')
hold on
plot(impacts(:,1),0,'rx')
ylim([-0.5 1])
legend('estimated impacts','real impacts')

% only looking at impact detections that have score > 1
multiple_detect = find(scores > 1);
figure;
plot(final_pred(multiple_detect)./Fs_pcb,0.5,'bx')
hold on
plot(impacts(:,1),0,'rx')
ylim([-0.5 1])
legend('estimated impacts','real impacts')

%% extract features and classify

detected_starts = final_pred(multiple_detect);
est_impactN = length(detected_starts);
est_arrive_idx_all = zeros(est_impactN,sensorN);
est_last_idx_all = zeros(est_impactN,sensorN);
est_width = zeros(est_impactN,sensorN);
est_peak_mag = zeros(est_impactN,sensorN);
est_energy = zeros(est_impactN,sensorN);
est_cwt_energy = zeros(est_impactN,sensorN);
est_cwt_peak = zeros(est_impactN,sensorN);
est_overlapping = zeros(est_impactN,1);
overlap_thresh = 1000; % # indexes that's considered an overlapping impact
window_width = 0.35; % in (s) of how long to set initial lastidx of window

for seg = 2:length(segments)
    start_impact_idx = segments(seg-1);
    last_impact_idx = segments(seg)-1;
    start_time_fsr = impacts(start_impact_idx:last_impact_idx,1);
    pcb_start_idx = round((start_time_fsr(1)-0.5)*Fs_pcb);
    pcb_last_idx = min([round((start_time_fsr(end)+0.5)*Fs_pcb), length(wien_pcbD(:,1))]);
    est_idx = find(detected_starts < pcb_last_idx & detected_starts > pcb_start_idx);
    est_start_times = detected_starts(est_idx)./Fs_pcb;
    for s = 1:sensorN
        seg_pcb = wien_pcbD(pcb_start_idx:pcb_last_idx,s);
        [wt,f] = cwt(seg_pcb,Fs_pcb);
        valid_f_idx = find(freq_higher & f > freq_lower);
        cwt_mag = abs(wt(valid_f_idx,:));
        sum_cwt = sum(cwt_mag,1);
        sum_smooth_cwt = movmean(sum_cwt, 1000);
        figure; 
        plot(sum_smooth_cwt)
        hold on
        for i = 1:length(est_start_times)
            global_idx = est_idx(1) + i - 1;
            fsr_start_idx = round(est_start_times(i)*Fs_pcb);
            starti = fsr_start_idx - 1000;
            if i ~= length(est_start_times) % not last element
                if est_start_times(i+1)*Fs_pcb - est_start_times(i)*Fs_pcb > overlap_thresh
                    lasti = min([round(est_start_times(i+1)*Fs_pcb), starti + round(window_width*Fs_pcb)]);
                else
                    lasti = starti + round(window_width*Fs_pcb);
                end
            else
                lasti = min([pcb_last_idx,starti + round(window_width*Fs_pcb)]);
            end
            cwt_window = sum_smooth_cwt(starti-pcb_start_idx:lasti-pcb_start_idx);
            % find arrival idx with aic pick
            aic_idx = aic_pick(cwt_window, 'to_peak');
            arrive_idx = aic_idx + starti;
            if arrive_idx < fsr_start_idx
                arrive_idx = fsr_start_idx;
            end
            % last_idx is the 1st idx that goes under thresh after peak
            cwt_window = sum_smooth_cwt(arrive_idx:lasti-pcb_start_idx);
            [~,maxidx] = max(cwt_window);
            under_thresh_idx = find(cwt_window(maxidx:end) < 0.00065);
            if isempty(under_thresh_idx)
                last_idx = lasti;
            else
                last_idx = under_thresh_idx(1) + aic_idx - 1 + maxidx - 1 + starti;
            end
            % store vars
            est_arrive_idx_all(global_idx,s) = arrive_idx;
            est_last_idx_all(global_idx,s) = last_idx;
            est_width(global_idx,s) = last_idx - arrive_idx;
            cwt_window = sum_smooth_cwt(arrive_idx-pcb_start_idx:last_idx-pcb_start_idx);
            est_cwt_peak(global_idx,s) = max(cwt_window);
            est_cwt_energy(global_idx,s) = sum(abs(cwt_window).^2);
            pcb_window = wien_pcbD(arrive_idx:last_idx,s);
            est_peak_mag(global_idx,s) = max(pcb_window);
            est_energy(global_idx,s) = sum(abs(pcb_window).^2);
            % overlap check
            if i > 1
                if est_last_idx_all(global_idx-1,s)-5 > arrive_idx
                    est_overlapping(global_idx) = 1;
                    est_overlapping(global_idx-1) = 1;
                end
            end
            % sanity check plot
            plot(arrive_idx-pcb_start_idx,0,'rx','MarkerSize',8)
            plot(last_idx-pcb_start_idx,0,'gx','MarkerSize',10)
        end
    end
end

X_test = [est_width,est_cwt_peak,est_cwt_energy,est_peak_mag,est_energy];
save(processedfilepath,'X_test','detected_starts','estimated_impacts',...
    'final_pred','multiple_detect','est_impactN','est_width','est_cwt_peak',...
    'est_cwt_energy','est_peak_mag','est_energy','est_overlapping','-append')


%% predict with tree, store everything in est_impacts

tree = fitctree(X_train,Y_train);
Y_predict = predict(tree,X_test);

est_impacts = zeros(length(Y_predict),10,6);
% sensors
for s = 1:6
    est_impacts(:,:,s) = [detected_starts,Y_predict,est_arrive_idx_all(:,s),est_last_idx_all(:,s),est_width(:,s),...
        est_cwt_peak(:,s),est_cwt_energy(:,s),est_peak_mag(:,s),est_energy(:,s),est_overlapping];
end

plot_footfall_labels(est_impacts(:,2,1),impacts,est_impacts(:,1,1),Fs_pcb)
labeling_success_rate(impacts, est_impacts(:,1,1), est_impacts(:,2,1), Fs_pcb, 0.1)

save(processedfilepath,'est_impacts','-append')

%% recursively iterate through to fix labeling using step time

% find segment end times
detected_seg_starts = [detected_starts(1)/Fs_pcb]; % start times of segments
for i = 2:length(detected_starts)
    if (detected_starts(i) - detected_starts(i-1))/Fs_pcb > 1.5 
        detected_seg_starts(end+1) = detected_starts(i)/Fs_pcb;
    end
end
detected_seg_starts(end+1) = detected_starts(end)/Fs_pcb + 0.5;

detected_seg_ends = []; % start times of segments
for i = 2:length(detected_starts)
    if (detected_starts(i) - detected_starts(i-1))/Fs_pcb > 1.5 
        detected_seg_ends(end+1) = detected_starts(i-1)/Fs_pcb;
    end
end
detected_seg_ends(end+1) = detected_starts(end)/Fs_pcb;

% first fix any big steps
% big o steps
o_idx = find(Y_predict == 1);
o_wrong = false;
for i = 2:length(o_idx)
    curr_time = detected_starts(o_idx(i))/Fs_pcb;
    prev_time = detected_starts(o_idx(i-1))/Fs_pcb;
    if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_time,1) % in same segment
        if curr_time - prev_time > 1.75*step_o1
            o_wrong = true;
            times_inbetween = detected_starts(o_idx(i-1):o_idx(i))./Fs_pcb - prev_time;
            [~,minidx] = min(abs(times_inbetween - step_o1));
            Y_predict(o_idx(i-1) + minidx - 1) = 1;
            o_idx = find(Y_predict == 1);
            break;
        end
    end
end
while o_wrong
    o_wrong = false;
    for i = 2:length(o_idx)
        curr_time = detected_starts(o_idx(i))/Fs_pcb;
        prev_time = detected_starts(o_idx(i-1))/Fs_pcb;
        if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_time,1) % in same segment
            if curr_time - prev_time > 1.75*step_o1
                o_wrong = true;
                times_inbetween = detected_starts(o_idx(i-1):o_idx(i))./Fs_pcb - prev_time;
                [~,minidx] = min(abs(times_inbetween - step_o1));
                Y_predict(o_idx(i-1) + minidx - 1) = 1;
                o_idx = find(Y_predict == 1);
                break;
            end
        end
    end
end

% big x steps
o_idx = find(Y_predict == 1);
x_idx = find(Y_predict == 2);
x_wrong = false;
for i = 2:length(x_idx)
    curr_time = detected_starts(x_idx(i))/Fs_pcb;
    prev_time = detected_starts(x_idx(i-1))/Fs_pcb;
    if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_time,1) % in same segment
        if curr_time - prev_time > 1.75*step_x1
            x_wrong = true;
            times_inbetween = detected_starts(x_idx(i-1):x_idx(i))./Fs_pcb - prev_time;
            [~,minidx] = min(abs(times_inbetween - step_x1));
            % check if this should be replaced with x or made overlapping
            curr_time
            minidx
            curr_idx = x_idx(i-1) + minidx - 1;
            curr_o_idx = find(o_idx == curr_idx);
            next_o_time = detected_starts(o_idx(curr_o_idx)+1)./Fs_pcb;
            prev_o_time = detected_starts(o_idx(curr_o_idx)-1)./Fs_pcb;
            curr_o_time = detected_starts(o_idx(curr_o_idx))./Fs_pcb;
            if abs((next_o_time-prev_o_time) - step_o1) < abs((next_o_time-curr_o_time) - step_o1)
                % replace 
                Y_predict(curr_idx) = 2;
                x_idx = find(Y_predict == 2);
            else
                % overlapping
                Y_predict = [Y_predict(1:curr_idx);2;Y_predict(curr_idx+1:end)];
                detected_starts = [detected_starts(1:curr_idx),detected_starts(curr_idx:end)];
            end
            x_idx = find(Y_predict == 2);
            break;
        end
    end
end
while x_wrong
    x_wrong = false;
    for i = 2:length(x_idx)
        curr_time = detected_starts(x_idx(i))/Fs_pcb;
        prev_time = detected_starts(x_idx(i-1))/Fs_pcb;
        if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_time,1) % in same segment
            if curr_time - prev_time > 1.75*step_x1
                x_wrong = true;
                times_inbetween = detected_starts(x_idx(i-1):x_idx(i))./Fs_pcb - prev_time;
                [~,minidx] = min(abs(times_inbetween - step_x1));
                % check if this should be replaced with x or made overlapping
                curr_time
                minidx
                curr_idx = x_idx(i-1) + minidx - 1;
                curr_o_idx = find(o_idx == curr_idx);
                next_o_time = detected_starts(o_idx(curr_o_idx)+1)./Fs_pcb;
                prev_o_time = detected_starts(o_idx(curr_o_idx)-1)./Fs_pcb;
                curr_o_time = detected_starts(o_idx(curr_o_idx))./Fs_pcb;
                if abs((next_o_time-prev_o_time) - step_o1) < abs((next_o_time-curr_o_time) - step_o1)
                    % replace 
                    Y_predict(curr_idx) = 2;
                    x_idx = find(Y_predict == 2);
                else
                    % overlapping
                    Y_predict = [Y_predict(1:curr_idx);2;Y_predict(curr_idx+1:end)];
                    detected_starts = [detected_starts(1:curr_idx),detected_starts(curr_idx:end)];
                end
                x_idx = find(Y_predict == 2);
                break;
            end
        end
    end
end

plot_footfall_labels(Y_predict,impacts,detected_starts,Fs_pcb,Fs_fsr)
labeling_success_rate(impacts, detected_starts, Y_predict, Fs_pcb, Fs_fsr, 0.1)

%% delete any small steps for o
small_scaler = 0.4;
% first fix any big steps
o_idx = find(Y_predict == 1);
x_idx = find(Y_predict == 2);
o_wrong = false;
for i = 2:length(o_idx)
    curr_time = detected_starts(o_idx(i))/Fs_pcb;
    prev_time = detected_starts(o_idx(i-1))/Fs_pcb;
    if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_time,1) % in same segment
        if curr_time - prev_time < small_scaler*step_o1
            o_wrong = true;
            % decide if remove or relabel for both impacts
            x_prev = find(detected_starts(x_idx) <= curr_time*Fs_pcb);
            x_prev = x_prev(end);
            x_prev_time = detected_starts(x_idx(x_prev))/Fs_pcb;
            x_next = find(detected_starts(x_idx) > curr_time*Fs_pcb);
            x_next = x_next(1);
            x_next_time = detected_starts(x_idx(x_next))/Fs_pcb;

            % check that prev and next x are in same seg
            if find(detected_seg_starts > x_next_time,1) == find(detected_seg_starts > curr_time,1) % in same segment
                if find(detected_seg_starts > x_prev_time,1) == find(detected_seg_starts > curr_time,1)
                    curr_or_prev = 0;
                    curr_margin = abs(x_next_time - curr_time - step_x1) + abs(curr_time - x_prev_time - step_x1);
                    prev_margin = abs(x_next_time - prev_time - step_x1) + abs(prev_time - x_prev_time - step_x1);
                    % check if better to replace
                    if abs(x_next_time - x_prev_time - step_x1) > abs(x_next_time - curr_time - step_x1)
                        if abs(x_next_time - x_prev_time - step_x1) > abs(curr_time - x_prev_time - step_x1)
                            curr_or_prev = 1;
                        end
                    end
                    if abs(x_next_time - x_prev_time - step_x1) > abs(x_next_time - prev_time - step_x1)
                        if abs(x_next_time - x_prev_time - step_x1) > abs(prev_time - x_prev_time - step_x1)
                            if curr_or_prev == 1 & prev_margin < curr_margin
                                curr_or_prev = 2;
                            elseif curr_or_prev == 0
                                curr_or_prev = 2;
                            end
                        end
                    end
                else % only next x is in same seg
                    curr_or_prev = 0;
                    replaced_curr_time = abs(x_next_time - curr_time - step_x1);
                    replaced_prev_time = abs(x_next_time - prev_time - step_x1);
                    if replaced_prev_time < replaced_curr_time
                        if replaced_prev_time < step_x1*0.3 
                            curr_or_prev = 2;
                        end
                    else
                        if replaced_curr_time < step_x1*0.3 
                            curr_or_prev = 1;
                        end
                    end
                end
            elseif find(detected_seg_starts > x_prev_time,1) == find(detected_seg_starts > curr_time,1)
                % only prev x is in same seg
                curr_or_prev = 0;
                replaced_curr_time = abs(curr_time - x_prev_time - step_x1);
                replaced_prev_time = abs(prev_time - x_prev_time - step_x1);
                if replaced_prev_time < replaced_curr_time
                    if replaced_prev_time < step_x1*0.3 
                        curr_or_prev = 2;
                    end
                else
                    if replaced_curr_time < step_x1*0.3 
                        curr_or_prev = 1;
                    end
                end
            end
            if curr_or_prev == 0 % remove impact
                next_o_time = detected_starts(o_idx(i+1))/Fs_pcb;
                prev_o_time = detected_starts(o_idx(i-2))/Fs_pcb;
                if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > next_o_time,1) % in same segment
                    if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_o_time,1) % in same segment
                        curr_o_margin = abs(next_o_time - curr_time - step_o1) + abs(curr_time - prev_o_time - step_o1);
                        prev_o_margin = abs(next_o_time - prev_time - step_o1) + abs(prev_time - prev_o_time - step_o1);

                    else
                        curr_o_margin = abs(next_o_time - curr_time - step_o1);
                        prev_o_margin = abs(next_o_time - prev_time - step_o1);
                    end
                elseif find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_o_time,1) % in same segment
                    curr_o_margin = abs(curr_time - prev_o_time - step_o1);
                    prev_o_margin = abs(prev_time - prev_o_time - step_o1);
                end
                if curr_o_margin > prev_o_margin % delete current impact
                    detected_starts(o_idx(i)) = [];
                    Y_predict(o_idx(i)) = [];
                else % delete prev impact
                    detected_starts(o_idx(i-1)) = [];
                    Y_predict(o_idx(i-1)) = [];
                end
            elseif curr_or_prev == 1 % replace current impact
                Y_predict(o_idx(i)) = 2;
            else % replace previous impact
                Y_predict(o_idx(i-1)) = 2;
            end
            o_idx = find(Y_predict == 1);
            x_idx = find(Y_predict == 2);
            break;
        end
    end
end
while o_wrong
    o_wrong = false;
    for i = 2:length(o_idx)
        curr_time = detected_starts(o_idx(i))/Fs_pcb;
        prev_time = detected_starts(o_idx(i-1))/Fs_pcb;
        if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_time,1) % in same segment
            if curr_time - prev_time < small_scaler*step_o1
                o_wrong = true;
                % decide if remove or relabel for both impacts
                x_prev = find(detected_starts(x_idx) <= curr_time*Fs_pcb);
                x_prev = x_prev(end);
                x_prev_time = detected_starts(x_idx(x_prev))/Fs_pcb;
                x_next = find(detected_starts(x_idx) > curr_time*Fs_pcb);
                x_next = x_next(1);
                x_next_time = detected_starts(x_idx(x_next))/Fs_pcb;

                % check that prev and next x are in same seg
                if find(detected_seg_starts > x_next_time,1) == find(detected_seg_starts > curr_time,1) % in same segment
                    if find(detected_seg_starts > x_prev_time,1) == find(detected_seg_starts > curr_time,1)
                        curr_or_prev = 0;
                        curr_margin = abs(x_next_time - curr_time - step_x1) + abs(curr_time - x_prev_time - step_x1);
                        prev_margin = abs(x_next_time - prev_time - step_x1) + abs(prev_time - x_prev_time - step_x1);
                        % check if better to replace
                        if abs(x_next_time - x_prev_time - step_x1) > abs(x_next_time - curr_time - step_x1)
                            if abs(x_next_time - x_prev_time - step_x1) > abs(curr_time - x_prev_time - step_x1)
                                curr_or_prev = 1;
                            end
                        end
                        if abs(x_next_time - x_prev_time - step_x1) > abs(x_next_time - prev_time - step_x1)
                            if abs(x_next_time - x_prev_time - step_x1) > abs(prev_time - x_prev_time - step_x1)
                                if curr_or_prev == 1 & prev_margin < curr_margin
                                    curr_or_prev = 2;
                                elseif curr_or_prev == 0
                                    curr_or_prev = 2;
                                end
                            end
                        end
                    else % only next x is in same seg
                        curr_or_prev = 0;
                        replaced_curr_time = abs(x_next_time - curr_time - step_x1);
                        replaced_prev_time = abs(x_next_time - prev_time - step_x1);
                        if replaced_prev_time < replaced_curr_time
                            if replaced_prev_time < step_x1*0.3 
                                curr_or_prev = 2;
                            end
                        else
                            if replaced_curr_time < step_x1*0.3 
                                curr_or_prev = 1;
                            end
                        end
                    end
                elseif find(detected_seg_starts > x_prev_time,1) == find(detected_seg_starts > curr_time,1)
                    % only prev x is in same seg
                    curr_or_prev = 0;
                    replaced_curr_time = abs(curr_time - x_prev_time - step_x1);
                    replaced_prev_time = abs(prev_time - x_prev_time - step_x1);
                    if replaced_prev_time < replaced_curr_time
                        if replaced_prev_time < step_x1*0.3 
                            curr_or_prev = 2;
                        end
                    else
                        if replaced_curr_time < step_x1*0.3 
                            curr_or_prev = 1;
                        end
                    end
                end
                if curr_or_prev == 0 % remove impact
                    next_o_time = detected_starts(o_idx(i+1))/Fs_pcb;
                    prev_o_time = detected_starts(o_idx(i-2))/Fs_pcb;
                    if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > next_o_time,1) % in same segment
                        if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_o_time,1) % in same segment
                            curr_o_margin = abs(next_o_time - curr_time - step_o1) + abs(curr_time - prev_o_time - step_o1);
                            prev_o_margin = abs(next_o_time - prev_time - step_o1) + abs(prev_time - prev_o_time - step_o1);
                            
                        else
                            curr_o_margin = abs(next_o_time - curr_time - step_o1);
                            prev_o_margin = abs(next_o_time - prev_time - step_o1);
                        end
                    elseif find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_o_time,1) % in same segment
                        curr_o_margin = abs(curr_time - prev_o_time - step_o1);
                        prev_o_margin = abs(prev_time - prev_o_time - step_o1);
                    end
                    if curr_o_margin > prev_o_margin % delete current impact
                        detected_starts(o_idx(i)) = [];
                        Y_predict(o_idx(i)) = [];
                    else % delete prev impact
                        detected_starts(o_idx(i-1)) = [];
                        Y_predict(o_idx(i-1)) = [];
                    end
                elseif curr_or_prev == 1 % replace current impact
                    Y_predict(o_idx(i)) = 2;
                else % replace previous impact
                    Y_predict(o_idx(i-1)) = 2;
                end
                o_idx = find(Y_predict == 1);
                x_idx = find(Y_predict == 2);
                break;
            end
        end
    end
end

% delete any small steps for x
small_scaler = 0.4;
% first fix any big steps
o_idx = find(Y_predict == 1);
x_idx = find(Y_predict == 2);
x_wrong = false;
for i = 2:length(x_idx)
    curr_time = detected_starts(x_idx(i))/Fs_pcb;
    prev_time = detected_starts(x_idx(i-1))/Fs_pcb;
    if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_time,1) % in same segment
        if curr_time - prev_time < small_scaler*step_x1
            x_wrong = true;
            % decide if remove or relabel for both impacts
            o_prev = find(detected_starts(o_idx) <= curr_time*Fs_pcb);
            if ~isempty(o_prev)
                o_prev = o_prev(end);
                o_prev_time = detected_starts(o_idx(o_prev))/Fs_pcb;
            else
                o_prev_time = detected_starts(end)/Fs_pcb + 1;
            end
            o_next = find(detected_starts(o_idx) > curr_time*Fs_pcb);
            o_next = o_next(1);
            o_next_time = detected_starts(o_idx(o_next))/Fs_pcb;

            % check that prev and next x are in same seg
            if find(detected_seg_starts > o_next_time,1) == find(detected_seg_starts > curr_time,1) % in same segment
                if find(detected_seg_starts > o_prev_time,1) == find(detected_seg_starts > curr_time,1)
                    curr_or_prev = 0;
                    curr_margin = abs(o_next_time - curr_time - step_o1) + abs(curr_time - o_prev_time - step_o1);
                    prev_margin = abs(o_next_time - prev_time - step_o1) + abs(prev_time - o_prev_time - step_o1);
                    % check if better to replace
                    if abs(o_next_time - o_prev_time - step_o1) > abs(o_next_time - curr_time - step_o1)
                        if abs(o_next_time - o_prev_time - step_o1) > abs(curr_time - o_prev_time - step_o1)
                            curr_or_prev = 1;
                        end
                    end
                    if abs(o_next_time - o_prev_time - step_o1) > abs(o_next_time - prev_time - step_o1)
                        if abs(o_next_time - o_prev_time - step_o1) > abs(prev_time - o_prev_time - step_o1)
                            if curr_or_prev == 1 & prev_margin < curr_margin
                                curr_or_prev = 2;
                            elseif curr_or_prev == 0
                                curr_or_prev = 2;
                            end
                        end
                    end
                else % only next x is in same seg
                    curr_or_prev = 0;
                    replaced_curr_time = abs(o_next_time - curr_time - step_o1);
                    replaced_prev_time = abs(o_next_time - prev_time - step_o1);
                    if replaced_prev_time < replaced_curr_time
                        if replaced_prev_time < step_o1*0.3 
                            curr_or_prev = 2;
                        end
                    else
                        if replaced_curr_time < step_o1*0.3 
                            curr_or_prev = 1;
                        end
                    end
                end
            elseif find(detected_seg_starts > o_prev_time,1) == find(detected_seg_starts > curr_time,1)
                % only prev x is in same seg
                curr_or_prev = 0;
                replaced_curr_time = abs(curr_time - o_prev_time - step_o1);
                replaced_prev_time = abs(prev_time - o_prev_time - step_o1);
                if replaced_prev_time < replaced_curr_time
                    if replaced_prev_time < step_o1*0.3 
                        curr_or_prev = 2;
                    end
                else
                    if replaced_curr_time < step_o1*0.3 
                        curr_or_prev = 1;
                    end
                end
            end
            if curr_or_prev == 0 % remove impact
                next_x_time = detected_starts(x_idx(i+1))/Fs_pcb;
                prev_x_time = detected_starts(x_idx(i-2))/Fs_pcb;
                if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > next_x_time,1) % in same segment
                    if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_x_time,1) % in same segment
                        curr_x_margin = abs(next_x_time - curr_time - step_x1) + abs(curr_time - prev_x_time - step_x1);
                        prev_x_margin = abs(next_x_time - prev_time - step_x1) + abs(prev_time - prev_x_time - step_x1);

                    else
                        curr_x_margin = abs(next_x_time - curr_time - step_x1);
                        prev_x_margin = abs(next_x_time - prev_time - step_x1);
                    end
                elseif find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_x_time,1) % in same segment
                    curr_x_margin = abs(curr_time - prev_x_time - step_x1);
                    prev_x_margin = abs(prev_time - prev_x_time - step_x1);
                end
                if curr_x_margin > prev_x_margin % delete current impact
                    detected_starts(x_idx(i)) = [];
                    Y_predict(x_idx(i)) = [];
                else % delete prev impact
                    detected_starts(x_idx(i-1)) = [];
                    Y_predict(x_idx(i-1)) = [];
                end
            elseif curr_or_prev == 1 % replace current impact
                Y_predict(x_idx(i)) = 1;
            else % replace previous impact
                Y_predict(x_idx(i-1)) = 1;
            end
            o_idx = find(Y_predict == 1);
            x_idx = find(Y_predict == 2);
            break;
        end
    end
end
while x_wrong
    x_wrong = false;
    for i = 2:length(x_idx)
        curr_time = detected_starts(x_idx(i))/Fs_pcb;
        prev_time = detected_starts(x_idx(i-1))/Fs_pcb;
        if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_time,1) % in same segment
            if curr_time - prev_time < small_scaler*step_x1
                x_wrong = true;
                % decide if remove or relabel for both impacts
                o_prev = find(detected_starts(o_idx) <= curr_time*Fs_pcb);
                o_prev = o_prev(end);
                o_prev_time = detected_starts(o_idx(o_prev))/Fs_pcb;
                o_next = find(detected_starts(o_idx) > curr_time*Fs_pcb);
                o_next = o_next(1);
                o_next_time = detected_starts(o_idx(o_next))/Fs_pcb;

                % check that prev and next x are in same seg
                if find(detected_seg_starts > o_next_time,1) == find(detected_seg_starts > curr_time,1) % in same segment
                    if find(detected_seg_starts > o_prev_time,1) == find(detected_seg_starts > curr_time,1)
                        curr_or_prev = 0;
                        curr_margin = abs(o_next_time - curr_time - step_o1) + abs(curr_time - o_prev_time - step_o1);
                        prev_margin = abs(o_next_time - prev_time - step_o1) + abs(prev_time - o_prev_time - step_o1);
                        % check if better to replace
                        if abs(o_next_time - o_prev_time - step_o1) > abs(o_next_time - curr_time - step_o1)
                            if abs(o_next_time - o_prev_time - step_o1) > abs(curr_time - o_prev_time - step_o1)
                                curr_or_prev = 1;
                            end
                        end
                        if abs(o_next_time - o_prev_time - step_o1) > abs(o_next_time - prev_time - step_o1)
                            if abs(o_next_time - o_prev_time - step_o1) > abs(prev_time - o_prev_time - step_o1)
                                if curr_or_prev == 1 & prev_margin < curr_margin
                                    curr_or_prev = 2;
                                elseif curr_or_prev == 0
                                    curr_or_prev = 2;
                                end
                            end
                        end
                    else % only next x is in same seg
                        curr_or_prev = 0;
                        replaced_curr_time = abs(o_next_time - curr_time - step_o1);
                        replaced_prev_time = abs(o_next_time - prev_time - step_o1);
                        if replaced_prev_time < replaced_curr_time
                            if replaced_prev_time < step_o1*0.3 
                                curr_or_prev = 2;
                            end
                        else
                            if replaced_curr_time < step_o1*0.3 
                                curr_or_prev = 1;
                            end
                        end
                    end
                elseif find(detected_seg_starts > o_prev_time,1) == find(detected_seg_starts > curr_time,1)
                    % only prev x is in same seg
                    curr_or_prev = 0;
                    replaced_curr_time = abs(curr_time - o_prev_time - step_o1);
                    replaced_prev_time = abs(prev_time - o_prev_time - step_o1);
                    if replaced_prev_time < replaced_curr_time
                        if replaced_prev_time < step_o1*0.3 
                            curr_or_prev = 2;
                        end
                    else
                        if replaced_curr_time < step_o1*0.3 
                            curr_or_prev = 1;
                        end
                    end
                end
                if curr_or_prev == 0 % remove impact
                    next_x_time = detected_starts(x_idx(i+1))/Fs_pcb;
                    prev_x_time = detected_starts(x_idx(i-2))/Fs_pcb;
                    if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > next_x_time,1) % in same segment
                        if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_x_time,1) % in same segment
                            curr_x_margin = abs(next_x_time - curr_time - step_x1) + abs(curr_time - prev_x_time - step_x1);
                            prev_x_margin = abs(next_x_time - prev_time - step_x1) + abs(prev_time - prev_x_time - step_x1);

                        else
                            curr_x_margin = abs(next_x_time - curr_time - step_x1);
                            prev_x_margin = abs(next_x_time - prev_time - step_x1);
                        end
                    elseif find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_x_time,1) % in same segment
                        curr_x_margin = abs(curr_time - prev_x_time - step_x1);
                        prev_x_margin = abs(prev_time - prev_x_time - step_x1);
                    end
                    if curr_x_margin > prev_x_margin % delete current impact
                        detected_starts(x_idx(i)) = [];
                        Y_predict(x_idx(i)) = [];
                    else % delete prev impact
                        detected_starts(x_idx(i-1)) = [];
                        Y_predict(x_idx(i-1)) = [];
                    end
                elseif curr_or_prev == 1 % replace current impact
                    Y_predict(x_idx(i)) = 1;
                else % replace previous impact
                    Y_predict(x_idx(i-1)) = 1;
                end
                o_idx = find(Y_predict == 1);
                x_idx = find(Y_predict == 2);
                break;
            end
        end
    end
end

plot_footfall_labels(Y_predict,impacts,detected_starts,Fs_pcb,Fs_fsr)
labeling_success_rate(impacts, detected_starts, Y_predict, Fs_pcb, Fs_fsr, 0.1)

%% finally, check edges of segments to see if any can be overlapping

overlap_wrong = true;
o_idx = find(Y_predict == 1);
x_idx = find(Y_predict == 2);
while overlap_wrong
    overlap_wrong = false;
    for i = 1:length(Y_predict)
        if ~isempty(find(round(detected_seg_starts.*Fs_pcb./10) == round(detected_starts(i)/10),1))
            % it's the start of the segment
            o_first_idx = find(detected_starts(o_idx) >= detected_starts(i),1);
            o_first = detected_starts(o_idx(o_first_idx));
            x_first_idx = find(detected_starts(x_idx) >= detected_starts(i),1);
            x_first = detected_starts(x_idx(x_first_idx));
            changed = false;
            if o_first < x_first
                % o happened first
                o_last_idx = find(detected_starts(o_idx) <= x_first & detected_starts(o_idx) >= o_first);
                % loop through all o's before x to see if overlapping x fit
                for j = 1:length(o_last_idx)
                    curr_idx = o_idx(o_last_idx(end-j + 1));
                    o_last = detected_starts(curr_idx);
                    if abs((x_first - o_last)/Fs_pcb - step_x1) < step_x1*0.36
                        % better to make overlapping
                        Y_predict = [Y_predict(1:curr_idx);2;Y_predict(curr_idx+1:end)];
                        detected_starts = [detected_starts(1:curr_idx),detected_starts(curr_idx:end)];
                        overlap_wrong = true;
                        changed = true;
                        o_idx = find(Y_predict == 1);
                        x_idx = find(Y_predict == 2);
                        break;
                    end
                end
            elseif o_first > x_first
                % x happened first
                x_last_idx = find(detected_starts(x_idx) <= o_first & detected_starts(x_idx) >= x_first);
                % loop through all o's before x to see if overlapping x fit
                for j = 1:length(x_last_idx)
                    curr_idx = x_idx(x_last_idx(end-j + 1));
                    x_last = detected_starts(curr_idx);
                    if abs((o_first - x_last)/Fs_pcb - step_o1) < step_o1*0.36
                        % better to make overlapping
                        Y_predict = [Y_predict(1:curr_idx);1;Y_predict(curr_idx+1:end)];
                        detected_starts = [detected_starts(1:curr_idx),detected_starts(curr_idx:end)];
                        overlap_wrong = true;
                        changed = true;
                        o_idx = find(Y_predict == 1);
                        x_idx = find(Y_predict == 2);
                        break;
                    end
                end
            end
            if changed
                break;
            end
        end
        % check for ends of segments
        if ~isempty(find(round(detected_seg_ends.*Fs_pcb./10) == round(detected_starts(i)/10),1))
            % it's the end of the segment
            o_last_idx = find(detected_starts(o_idx) <= detected_starts(i));
            o_last = detected_starts(o_idx(o_last_idx(end)));
            x_last_idx = find(detected_starts(x_idx) <= detected_starts(i));
            x_last = detected_starts(x_idx(x_last_idx(end)));
            changed = false;
            if o_last > x_last
                % o happened last
                o_last_idx = find(detected_starts(o_idx) >= x_last & detected_starts(o_idx) <= o_last);
                % loop through all o's before x to see if overlapping x fit
                for j = 1:length(o_last_idx)
                    curr_idx = o_idx(o_last_idx(j));
                    curr_last = detected_starts(curr_idx);
                    if abs((curr_last - x_last)/Fs_pcb - step_x1) < step_x1*0.36
                        % better to make overlapping
                        Y_predict = [Y_predict(1:curr_idx);2;Y_predict(curr_idx+1:end)];
                        detected_starts = [detected_starts(1:curr_idx),detected_starts(curr_idx:end)];
                        overlap_wrong = true;
                        changed = true;
                        o_idx = find(Y_predict == 1);
                        x_idx = find(Y_predict == 2);
                        break;
                    end
                end
            elseif o_last < x_last
                % x happened last
                x_last_idx = find(detected_starts(x_idx) >= o_last & detected_starts(x_idx) <= x_last);
                % loop through all o's before x to see if overlapping x fit
                for j = 1:length(x_last_idx)
                    curr_idx = x_idx(x_last_idx(j));
                    curr_last = detected_starts(curr_idx);
                    if abs((curr_last - o_last)/Fs_pcb - step_o1) < step_o1*0.36
                        % better to make overlapping
                        Y_predict = [Y_predict(1:curr_idx);1;Y_predict(curr_idx+1:end)];
                        detected_starts = [detected_starts(1:curr_idx),detected_starts(curr_idx:end)];
                        overlap_wrong = true;
                        changed = true;
                        o_idx = find(Y_predict == 1);
                        x_idx = find(Y_predict == 2);
                        break;
                    end
                end
            end
            if changed
                break;
            end
        end
    end
end

plot_footfall_labels(Y_predict,impacts,detected_starts,Fs_pcb,Fs_fsr)
labeling_success_rate(impacts, detected_starts, Y_predict, Fs_pcb, Fs_fsr, 0.1)

% plots for manually removing for next section
figure;
est_o = find(Y_predict == 1);
est_x = find(Y_predict == 2);
plot(detected_starts(est_o),est_o,'ro')
hold on
plot(detected_starts(est_x),est_x,'rx')

figure;
real_o = find(impacts(:,4) == 1 | impacts(:,4) == 2);
real_x = find(impacts(:,4) == 3 | impacts(:,4) == 4);
plot(impacts(real_o,1)./Fs_fsr,real_o,'bo')
hold on
plot(impacts(real_x,1)./Fs_fsr,real_x,'bx')



%% manually label parts in start/end of segments to find ultimate success rate
remove_est_idx = [];
remove_real_idx = [];

clean_detected_starts = detected_starts;
clean_Y_predict = Y_predict;
clean_impacts = impacts;

clean_detected_starts(remove_est_idx) = [];
clean_Y_predict(remove_est_idx) =[];
clean_impacts(remove_real_idx,:) = [];

plot_footfall_labels(clean_Y_predict,clean_impacts,clean_detected_starts,Fs_pcb,Fs_fsr)
labeling_success_rate(clean_impacts, clean_detected_starts, clean_Y_predict, Fs_pcb, Fs_fsr, 0.1)


%% save vars
save(processedfilepath,'detected_seg_starts','detected_seg_ends','Y_predict',...
    'detected_starts','clean_detected_starts','clean_impacts','clean_Y_predict','-append')
