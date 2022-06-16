%% use decision tree to classify impacts
% load training data
load('C:\Users\katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Jenny 1\ProcessedData\both_limp1.mat')
tree = fitctree(X_train,Y_train);

% load testing data
processedfilepath = 'C:\Users\katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Jenny 1\ProcessedData\both_limp2.mat';
load(processedfilepath)

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

% only looking at impact detections that have score > 1
multiple_detect = find(scores > 1);
figure;
plot(final_pred(multiple_detect)./Fs_pcb,0.5,'bx')
hold on
plot(impacts(:,1),0,'rx')
ylim([-0.5 1])

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
% save(processedfilepath,'X_test','detected_starts','estimated_impacts',...
%     'final_pred','multiple_detect','est_impactN','est_width','est_cwt_peak',...
%     'est_cwt_energy','est_peak_mag','est_energy','est_overlapping','-append')


%% predict with tree, store everything in est_impacts

Y_predict = predict(tree,X_test);

est_impacts = [detected_starts.',Y_predict,est_arrive_idx_all,est_last_idx_all,est_width,...
        est_cwt_peak,est_cwt_energy,est_peak_mag,est_energy,est_overlapping];

plot_footfall_labels(est_impacts(:,2),impacts,est_impacts(:,1),Fs_pcb)
title('Just classifier')
labeling_success_rate(impacts, est_impacts(:,1), est_impacts(:,2), Fs_pcb, 0.1)
labeling_success_rate(impacts, est_impacts(:,1), est_impacts(:,2), Fs_pcb, 0.05)

save(processedfilepath,'est_impacts','-append')

%% recursively iterate through to fix labeling using step time
% this is where it starts getting different from _processing3.m
% ######################################################################

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
old_est_impacts = est_impacts;

is_startBig = false;
startBig_est_impacts = est_impacts;
startSmall_est_impacts = est_impacts;
% big o steps
fixed = true;
while fixed
    [startBig_est_impacts,fixed] = big_step_labeling_o_limp(startBig_est_impacts,detected_seg_starts,step_o1,step_o2,Fs_pcb,true);
end
fixed = true;
while fixed
    [startSmall_est_impacts,fixed] = big_step_labeling_o_limp(startSmall_est_impacts,detected_seg_starts,step_o1,step_o2,Fs_pcb,false);
end
big_success_rate = labeling_success_rate(impacts, startBig_est_impacts(:,1), startBig_est_impacts(:,2), Fs_pcb, 0.05);
small_success_rate = labeling_success_rate(impacts, startSmall_est_impacts(:,1), startSmall_est_impacts(:,2), Fs_pcb, 0.05);
if big_success_rate > small_success_rate
    is_startBig = true;
    est_impacts = startBig_est_impacts;
else
    est_impacts = startSmall_est_impacts;
end

%% big x steps
fixed = true;
while fixed
    [est_impacts,fixed] = big_step_labeling_x_limp(est_impacts,detected_seg_starts,step_x1,step_o1,step_o2,is_startBig,Fs_pcb);
end

plot_footfall_labels(est_impacts(:,2),impacts,est_impacts(:,1),Fs_pcb)
title('After big step labeling')
labeling_success_rate(impacts, est_impacts(:,1), est_impacts(:,2), Fs_pcb, 0.1)
labeling_success_rate(impacts, est_impacts(:,1), est_impacts(:,2), Fs_pcb, 0.05)

%% delete any small steps
% small steps for o
fixed = true;
while fixed
    [est_impacts,fixed] = small_step_labeling_limp(est_impacts,detected_seg_starts,step_o1,step_o2,step_x1,1,2,is_startBig,Fs_pcb);
end

% small steps for x
fixed = true;
while fixed
    [est_impacts,fixed] = small_step_labeling_limp(est_impacts,detected_seg_starts,step_o1,step_o2,step_x1,2,1,is_startBig,Fs_pcb);
end

plot_footfall_labels(est_impacts(:,2),impacts,est_impacts(:,1),Fs_pcb)
title('After small step labeling')
labeling_success_rate(impacts, est_impacts(:,1), est_impacts(:,2), Fs_pcb, 0.1)
labeling_success_rate(impacts, est_impacts(:,1), est_impacts(:,2), Fs_pcb, 0.05)

%% finally, check edges of segments to see if any can be overlapping
% [est_impacts,~] = overlapping_step_labeling(est_impacts,detected_seg_starts,detected_seg_ends,step_o1,step_x1,Fs_pcb);
[est_impacts,~] = overlapping_step_labeling_limp(est_impacts,detected_seg_starts,detected_seg_ends,step_o1,step_o2,step_x1,is_startBig,Fs_pcb);

plot_footfall_labels(est_impacts(:,2),impacts,est_impacts(:,1),Fs_pcb)
title('After overlap labeling')
labeling_success_rate(impacts, est_impacts(:,1), est_impacts(:,2), Fs_pcb, 0.1)
labeling_success_rate(impacts, est_impacts(:,1), est_impacts(:,2), Fs_pcb, 0.05)

% %% repeat until no more fixes need to be made
% 
% count = 1;
% fixed = true;
% while fixed
%     [est_impacts,fixed] = big_step_labeling_o_limp(startBig_est_impacts,detected_seg_starts,step_o1,step_o2,Fs_pcb,is_startBig);
%     if fixed == true
%         count = count + 1;
%     end
% end
% fixed = true;
% while fixed
%     [est_impacts,fixed] = big_step_labeling_x_limp(est_impacts,detected_seg_starts,step_x1,step_o1,step_o2,is_startBig,Fs_pcb);
%     if fixed == true
%         count = count + 1;
%     end
% end
% fixed = true;
% while fixed
%     [est_impacts,fixed] = small_step_labeling_limp(est_impacts,detected_seg_starts,step_o1,step_o2,step_x1,1,2,is_startBig,Fs_pcb);
%     if fixed == true
%         count = count + 1;
%     end
% end
% fixed = true;
% while fixed
%     [est_impacts,fixed] = small_step_labeling_limp(est_impacts,detected_seg_starts,step_o1,step_o2,step_x1,2,1,is_startBig,Fs_pcb);
%     if fixed == true
%         count = count + 1;
%     end
% end
% [est_impacts,iter] = overlapping_step_labeling_limp(est_impacts,detected_seg_starts,detected_seg_ends,step_o1,step_o2,step_x1,is_startBig,Fs_pcb);
% %%
% while count ~= 1 | iter ~= 0
%     count = 1;
%     fixed = true;
%     while fixed
%         [est_impacts,fixed] = big_step_labeling_o_limp(startBig_est_impacts,detected_seg_starts,step_o1,step_o2,Fs_pcb,is_startBig);
%         if fixed == true
%             count = count + 1;
%         end
%     end
%     fixed = true;
%     while fixed
%         [est_impacts,fixed] = big_step_labeling_x_limp(est_impacts,detected_seg_starts,step_x1,step_o1,step_o2,is_startBig,Fs_pcb);
%         if fixed == true
%             count = count + 1;
%         end
%     end
%     fixed = true;
%     while fixed
%         [est_impacts,fixed] = small_step_labeling_limp(est_impacts,detected_seg_starts,step_o1,step_o2,step_x1,1,2,is_startBig,Fs_pcb);
%         if fixed == true
%             count = count + 1;
%         end
%     end
%     fixed = true;
%     while fixed
%         [est_impacts,fixed] = small_step_labeling_limp(est_impacts,detected_seg_starts,step_o1,step_o2,step_x1,2,1,is_startBig,Fs_pcb);
%         if fixed == true
%             count = count + 1;
%         end
%     end
%     [est_impacts,iter] = overlapping_step_labeling_limp(est_impacts,detected_seg_starts,detected_seg_ends,step_o1,step_o2,step_x1,is_startBig,Fs_pcb);
%     count
%     iter
% end
% 
% plot_footfall_labels(est_impacts(:,2),impacts,est_impacts(:,1),Fs_pcb)
% labeling_success_rate(impacts, est_impacts(:,1), est_impacts(:,2), Fs_pcb, 0.1)
% labeling_success_rate(impacts, est_impacts(:,1), est_impacts(:,2), Fs_pcb, 0.05)

%% plots for manually removing for next section
figure;
est_o = find(est_impacts(:,2) == 1);
est_x = find(est_impacts(:,2) == 2);
plot(est_impacts(est_o,1)./Fs_pcb,est_o,'ro')
hold on
plot(est_impacts(est_x,1)./Fs_pcb,est_x,'rx')

figure;
real_o = find(floor(impacts(:,2)/10) == 1);
real_x = find(floor(impacts(:,2)/10) == 2);
plot(impacts(real_o,1),real_o,'bo')
hold on
plot(impacts(real_x,1),real_x,'bx')

%% manually remove parts in start/end of segments to find ultimate success rate

remove_est_idx = [1,16,18,19,41,65,83,128,150,213,232,296,298,317,318]; % red graph
remove_real_idx = [16,17,80,82,108,234,236,259,275,277,279,281,326,324]; % blue graph

clean_est_impacts = est_impacts;
clean_impacts = impacts;

clean_est_impacts(remove_est_idx,:) = [];
clean_impacts(remove_real_idx,:) = [];

plot_footfall_labels(clean_est_impacts(:,2),clean_impacts,clean_est_impacts(:,1),Fs_pcb)
labeling_success_rate(clean_impacts, clean_est_impacts(:,1), clean_est_impacts(:,2), Fs_pcb, 0.1)
labeling_success_rate(clean_impacts, clean_est_impacts(:,1), clean_est_impacts(:,2), Fs_pcb, 0.05)

%% plots for manually removing for next section
figure;
est_o = find(clean_est_impacts(:,2) == 1);
est_x = find(clean_est_impacts(:,2) == 2);
plot(clean_est_impacts(est_o,1)./Fs_pcb,est_o,'ro')
hold on
plot(clean_est_impacts(est_x,1)./Fs_pcb,est_x,'rx')


%% manually add missing impacts
add_idx = [19,236];
add_label = [1,2];
replace_idx = [76,182,183,221];
replace_label = [1,2,1,1];

for i = 1:length(replace_idx)
    clean_est_impacts(replace_idx(i),2) = replace_label(i);
end
for i = 1:length(add_idx)
    curridx = add_idx(i) + i - 1; % to account for matrix increasing with more add
    new_array = clean_est_impacts(curridx,:);
    new_array(2) = add_label(i);
    clean_est_impacts = [clean_est_impacts(1:curridx,:);new_array;clean_est_impacts(curridx+1:end,:)];
end

plot_footfall_labels(clean_est_impacts(:,2),clean_impacts,clean_est_impacts(:,1),Fs_pcb)
labeling_success_rate(clean_impacts, clean_est_impacts(:,1), clean_est_impacts(:,2), Fs_pcb, 0.1)
labeling_success_rate(clean_impacts, clean_est_impacts(:,1), clean_est_impacts(:,2), Fs_pcb, 0.05)

%% delete section
est_seg_start = 90.9;
est_seg_end = 97.4;
seg_start = 91;
seg_end = 97.5;
delete_idx = find(clean_est_impacts(:,1)./Fs_pcb > est_seg_start & clean_est_impacts(:,1)./Fs_pcb < est_seg_end);
clean_est_impacts(delete_idx,:) = [];
delete_idx = find(clean_impacts(:,1) > seg_start & clean_impacts(:,1) < seg_end);
clean_impacts(delete_idx,:) = [];
plot_footfall_labels(clean_est_impacts(:,2),clean_impacts,clean_est_impacts(:,1),Fs_pcb)
labeling_success_rate(clean_impacts, clean_est_impacts(:,1), clean_est_impacts(:,2), Fs_pcb, 0.1)
labeling_success_rate(clean_impacts, clean_est_impacts(:,1), clean_est_impacts(:,2), Fs_pcb, 0.05)


%% save vars
save(processedfilepath,'detected_seg_starts','detected_seg_ends','est_impacts',...
    'clean_est_impacts','clean_impacts','-append')
disp(append("Saved as ", processedfilepath))
