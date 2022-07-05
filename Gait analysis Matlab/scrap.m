%% finding overlapping impact example for SFCWT

startidx = 840321;
lastidx = 917803;

for i = 1:6
    sensornum = i;

    clip = wien_pcbD(startidx:lastidx,sensornum);
    [wt,f] = cwt(clip,12800);
    valid_f_idx = find(f < freq_higher & f > freq_lower);
    cwt_freq = f(valid_f_idx);
    cwt_mag = abs(wt(valid_f_idx,:));
    sum_cwt = sum(cwt_mag,1);
    sum_smooth_cwt = movmean(sum_cwt, 1000);
    figure;
    subplot(2,1,1)
    plot(sum_smooth_cwt)
    hold on
    idx = find(impacts(:,1) > (startidx/12800) & impacts(:,1) < (lastidx/12800));
    plot(impacts(idx,1)*12800 - startidx,0,'rx','MarkerSize',10)

    subplot(2,1,2)
    pcbclip = pcbData(startidx:lastidx,sensornum);
    plot(pcbclip)
    hold on
    plot(impacts(idx,1)*12800 - startidx,0,'rx','MarkerSize',10)
    title(['start ' , num2str(startidx/12800) , ' sensor ' , num2str(sensornum)])
end


%% testing GMM derivation with dummy data

r1 = normrnd(3,1,[1,500]);
r2 = normrnd(4,2,[1,500]);

filename = 'C:\Users\katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\GMM_testing.csv';
writematrix([r1.';r2.'],filename)
% 6/17/22: it works! yay!

%% python GMM

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
        filename = 'C:\Users\katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\both_limp2_p1_GMM.csv';
        writematrix(step_times.',filename)
    else
        filename = 'C:\Users\katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\both_limp2_p2_GMM.csv';
        writematrix(step_times.',filename)
    end
end

%%
o_idx = find(est_impacts(:,2) == 1);
x_idx = find(est_impacts(:,2) == 2);
fixed = false;
for i = 6:length(x_idx)
    curr_time = est_impacts(x_idx(i),1)/Fs_pcb;
    prev_time = est_impacts(x_idx(i-1),1)/Fs_pcb;
    if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_time,1) % in same segment
        if curr_time - prev_time > 1.75*step_x1
            fixed = true;
            times_inbetween = est_impacts(x_idx(i-1):x_idx(i),1)./Fs_pcb - prev_time;
            [~,minidx] = min(abs(times_inbetween - step_x1));
            % check if this should be replaced with x or made overlapping
            curr_idx = x_idx(i-1) + minidx - 1;
            curr_o_idx = find(o_idx == curr_idx);
            next_o_time = est_impacts(o_idx(curr_o_idx)+1,1)./Fs_pcb;
            prev_o_time = est_impacts(o_idx(curr_o_idx)-1,1)./Fs_pcb;
            curr_o_time = est_impacts(o_idx(curr_o_idx),1)./Fs_pcb;
            if abs((next_o_time-prev_o_time) - step_o1) < abs((next_o_time-curr_o_time) - step_o1)
                % replace
                est_impacts(curr_idx,2) = 2;
            else
                % overlapping
                new_array = [est_impacts(curr_idx,1),2,est_impacts(curr_idx,3:end)];
                est_impacts = [est_impacts(1:curr_idx,:);new_array;est_impacts(curr_idx+1:end,:)];
            end
            break;
        end
    end
end

%% recursively iterate through to fix labeling using step time

detected_seg_ends = []; % start times of segments
for i = 2:length(detected_starts)
    if (detected_starts(i) - detected_starts(i-1))/Fs_pcb > 1.5 
        detected_seg_ends(end+1) = detected_starts(i-1)/Fs_pcb;
    end
end
detected_seg_ends(end+1) = detected_starts(end)/Fs_pcb;

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

                
                
                
            
    



%% use localization to check if proposed labels are correct 6/6/22

% this uses workspace data from below sections
w = 4;
curr_estimates = final_estimates_labels(curr_est,:);
impacttimes = impacts(segments(w-1):segments(w)-1,1)/Fs_fsr;

est_o = find(curr_estimates(:,2) == 1 | curr_estimates(:,2) == 2);
est_x = find(curr_estimates(:,2) == 3 | curr_estimates(:,2) == 4);

window_thresh = 0.0825; % seconds of how long an impact lasts
start_idx_thresh = 1000; % count 1000 idx before curr_time to do aicpick
freq_lower = 100; % Hz; frequency limits for cwt
freq_higher = 350;
estL = length(curr_estimates(:,1));

% initialize feature arrays
est_arrival_idx = zeros(estL,sensorN);
est_last_idx = zeros(estL,sensorN);
est_peak_idx = zeros(estL,sensorN);
est_peak_mag = zeros(estL,sensorN);
est_energy = zeros(estL,sensorN);
est_cwt_energy = zeros(estL,sensorN);
est_cwt_peak = zeros(estL,sensorN);
est_overlapping = zeros(estL,1); % value = 1 if two impacts overlap

for s = 1:sensorN
    % extract cwt of segment
    not_overlap = [];
    first_idx = round((curr_estimates(1,1)-0.5)*Fs_pcb);
    last_idx = min([round((curr_estimates(estL,1)+0.5)*Fs_pcb),length(wien_pcbD(:,1))]);
    seg_pcb = wien_pcbD(first_idx:last_idx,s);
    [wt,f] = cwt(seg_pcb,Fs_pcb);
    valid_f_idx = find(freq_higher & f > freq_lower);
    cwt_mag = abs(wt(valid_f_idx,:));
    sum_cwt = sum(cwt_mag,1);
    sum_smooth_cwt = movmean(sum_cwt, 800);
    pcb2cwt_idx = first_idx; % add this to result to get pcb idx
    % extract features from cwt sum
    for i = 1:estL
        % check for overlapping
        if i~= estL & round(curr_estimates(i,1),3) == round(curr_estimates(i+1,1),3)
            est_overlapping(i) = 1;
            est_overlapping(i+1) = 1;
        elseif est_overlapping(i) == 0
            not_overlap(end+1) = i;
            arrival_idx = round(curr_estimates(i,1)*Fs_pcb);
            if i ~= estL
                [lasti,nextimpact] = min([round(curr_estimates(i+1,1)*Fs_pcb),arrival_idx + round(0.46*Fs_pcb)]);
            else
                lasti = min([last_idx,arrival_idx + round(0.31*Fs_pcb)]);
                nextimpact = 2;
            end
            window = sum_smooth_cwt(arrival_idx-pcb2cwt_idx+1:lasti-pcb2cwt_idx+1);
            deriv_window = movmean(diff(window),1200);
            [~,minidx] = min(deriv_window);
            [~,maxidx] = max(deriv_window);
            est_arrival_idx(i,s) = arrival_idx;
            % window end time
            [~,zeroidx] = min(abs(deriv_window(minidx:end)));
            abs_last_idx = zeroidx + minidx - 1 + arrival_idx;
            est_last_idx(i,s) = abs_last_idx;
            % cwt energy
            if (abs_last_idx - arrival_idx) > 10 % sometimes peak is so small it's a straight line, then last idx is right after arrival idx
                window2 = wien_pcbD(arrival_idx:abs_last_idx,s);
                [wt2,f2] = cwt(window2,Fs_pcb);
                valid_f_idx2 = find(freq_higher & f2 > freq_lower);
                cwt_mag2 = abs(wt2(valid_f_idx2,:));
                sum_cwt2 = sum(cwt_mag2,1);
                sum_smooth_cwt2 = movmean(sum_cwt2, 800);
                est_cwt_energy(i,s) = sum(sum_smooth_cwt2.^2);
                % peak mag
                est_peak_mag(i,s) = max(abs(window2));
                % peak cwt mag
                est_cwt_peak(i,s) = max(sum_smooth_cwt2);
                % energy
                est_energy(i,s) = sum(abs(window2).^2);
            end
        end
    end
end

% sanity check
figure;
plot(seg_pcb)
hold on
plot(est_arrival_idx(:,6)-first_idx,0,'rx','MarkerSize',10)
plot(est_last_idx(:,6)-first_idx,0,'gx','MarkerSize',8)

% write matrix for localization
featureV = zeros(1,25);
for i = 1:estL
    if est_overlapping(i) == 0
        if curr_estimates(i,2) == 1
            curr_label = 1;
        elseif curr_estimates(i,2) == 2
            curr_label = 1;
        else
            curr_label = 2;
        end
        
        curr_feature = [curr_label,est_cwt_peak(i,:),est_cwt_energy(i,:),est_peak_mag(i,:),est_energy(i,:)];
        featureV(end+1,:) = curr_feature;
    end
end
featureV(1,:) = [];
writematrix(featureV,'C:\Users\katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\ExcelData\both_weight2_testing_060622.csv')

%% using the localization result to relabel 6/6/22

pythonfilename = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\ExcelData\both_weight2_testing_060622_results.csv';
T = readtable(pythonfilename);
A = table2array(T);

person1_loc = A(:,1);
person2_loc = A(:,2);
first_label = featureV(:,1);

person1_idx = find(first_label == 1);
person2_idx = find(first_label == 2);

new_curr_est = curr_estimates(find(est_overlapping == 0),:);
est_o = find(new_curr_est(:,2) == 1 | new_curr_est(:,2) == 2);
est_x = find(new_curr_est(:,2) == 3 | new_curr_est(:,2) == 4);

figure;
plot(person1_idx, person1_loc(person1_idx),'-o')
title('initial person 1 loc')
ylim([-4 4])
figure;
plot(person2_idx, person2_loc(person2_idx),'-o')
title('initial person 2 loc')
ylim([-4 4])

figure;
plot((impacttimes(o_impactidx)),0,'bo','MarkerSize',8)
hold on
plot((impacttimes(x_impactidx)),0,'bx','MarkerSize',8)
hold on
plot(new_curr_est(est_o,1),0.5,'ro')
plot(new_curr_est(est_x,1),0.5,'rx')
ylim([-0.5 1])

% correct localization
figure;
real_labels = [1,2,1,1,2,1,2,1,2,1,2,1,2,1,1,1,2];
new_person1_idx = find(real_labels == 1);
new_person2_idx = find(real_labels == 2);
plot(new_person1_idx, person1_loc(new_person1_idx),'-o')
title('correct person 1 loc')
ylim([-4 4])
figure;
plot(new_person2_idx, person2_loc(new_person2_idx),'-o')
title('correct person 2 loc')
ylim([-4 4])

% calculating all combinations
step_times = new_curr_est(:,1);
powerN = length(step_times);
low_powerN = round(powerN/3); % know they at least have to be 33% of total footsteps
high_powerN = round(2*powerN/3);
comb_len = 0;
for i = low_powerN:high_powerN
    comb_len = comb_len + nchoosek(powerN,i);
end
allcombinations = ones(comb_len, powerN);
scores = ones(comb_len, 1);
% get all unique perms
assign_idx = 1;
for i = low_powerN:high_powerN 
    curr_arr = ones(1,powerN);
    curr_arr(1:i) = 2;
    P = uniqueperms(curr_arr);
    Prows = size(P,1);
    allcombinations(assign_idx:assign_idx+Prows-1,:) = P;
    assign_idx = assign_idx + Prows;
end

for i = 1:comb_len
    curr_o = find(allcombinations(i,:) == 1);
    curr_x = find(allcombinations(i,:) == 2);
    locs_o = person1_loc(curr_o);
    locs_x = person2_loc(curr_x);
    scores(i) = sum(abs(diff(diff(locs_o)))) + sum(abs(diff(diff(locs_x))));
end

[~,minidx] = min(scores)
final_comb = allcombinations(minidx,:);
figure;
plot((impacttimes(o_impactidx)),0,'bo','MarkerSize',8)
hold on
plot((impacttimes(x_impactidx)),0,'bx','MarkerSize',8)
est_o = find(final_comb == 1);
est_x = find(final_comb == 2);
plot(step_times(est_o,1),0.5,'ro')
plot(step_times(est_x,1),0.5,'rx')
ylim([-0.5 1])




%% visualize to find an example where similar step times have one off labeling 6/3/22
w = 4;
start_time = impacts(segments(w-1),1)/Fs_fsr-0.15;
stop_time = impacts(segments(w)-1,1)/Fs_fsr+0.15;

start_idx_pcb = findTindex(start_time,pcbTime);
stop_idx_pcb = findTindex(stop_time,pcbTime);

impacttimes = impacts(segments(w-1):segments(w)-1,1)/Fs_fsr;


figure;
o_impactidx = find(impacts(segments(w-1):segments(w)-1,4) == 1 | impacts(segments(w-1):segments(w)-1,4) == 2);
x_impactidx = find(impacts(segments(w-1):segments(w)-1,4) == 3 | impacts(segments(w-1):segments(w)-1,4) == 4);
plot((impacttimes(o_impactidx)),0,'bo','MarkerSize',8)
hold on
plot((impacttimes(x_impactidx)),0,'bx','MarkerSize',8)

curr_est = find(final_estimates_labels(:,3) == w);
est_o = find(final_estimates_labels(curr_est,2) == 1 | final_estimates_labels(curr_est,2) == 2);
est_x = find(final_estimates_labels(curr_est,2) == 3 | final_estimates_labels(curr_est,2) == 4);
plot(final_estimates_labels(curr_est(est_o),1),0.5,'ro','MarkerSize',8)
plot(final_estimates_labels(curr_est(est_x),1),0.5,'rx','MarkerSize',8)
ylim([-0.5 1])


%%
load('C:\Users\katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\both_weight2.mat')
% figure;
% est_o = find(clean_final_estimates_labels(:,2) == 1);
% est_x = find(clean_final_estimates_labels(:,2) == 2);
% real_o = find(clean_ID_labels == 1);
% real_x = find(clean_ID_labels == 2);
% plot(clean_final_estimates_labels(est_o,1),0.2,'ro')
% hold on
% plot(clean_final_estimates_labels(est_x,1),0.2,'rx')
% plot(clean_impacts(real_o,1)./Fs_fsr,0,'bo')
% plot(clean_impacts(real_x,1)./Fs_fsr,0,'bx')
% ylim([-1,1])

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

%% cont from above
w = 4;
pkprom = [0.0015,0.0015,0.0015,0.0015,0.0015,0.0015];
pkheight = [0.006,0.016,0.0014,0.004,0.0127,0.009];
freq_lower = 100; % Hz; frequency limits for cwt
freq_higher = 350;
arrival_window = 1500; % 5000/Fs_pcb width of window to find arrival time
impact_time_thresh = 0.05; % if two impacts are 0.05s within each other, combine
correct_thresh = 0.05; % if predicted is within this thresh to real impact time, it's correct

start_time = impacts(segments(w-1),1)/Fs_fsr-0.15;
stop_time = impacts(segments(w)-1,1)/Fs_fsr+0.15;

start_idx_pcb = findTindex(start_time,pcbTime);
stop_idx_pcb = findTindex(stop_time,pcbTime);

impacttimes = impacts(segments(w-1):segments(w)-1,1)/Fs_fsr;

estimated_impacts = [0,0]; % [impact time, sensor N], for each segment

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
    impacttimes = impacts(segments(w-1):segments(w)-1,1)/Fs_fsr;
    o_impactidx = find(impacts(segments(w-1):segments(w)-1,4) == 1 | impacts(segments(w-1):segments(w)-1,4) == 2);
    x_impactidx = find(impacts(segments(w-1):segments(w)-1,4) == 3 | impacts(segments(w-1):segments(w)-1,4) == 4);
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
%             [minval, minidx] = min(sum_smooth_cwt(locs(i-1):locs(i)));
%             starti = max(locs(i) - arrival_window, locs(i-1)+minidx - arrival_window/5);
        end
        window = sum_smooth_cwt(starti:locs(i));
        arrive_idx = aic_pick(window, 'to_peak')+starti - 1;
%         if i ~= 1 & minval > 0.0001
%             arrive_idx = starti;
%         else
%             window = sum_smooth_cwt(starti:locs(i));
%             arrive_idx = aic_pick(window, 'to_peak')+starti;
%         end
        add_array = [arrive_idx, s];
        estimated_impacts = [estimated_impacts; add_array];
        plot(arrive_idx,0,'bx') 
    end
%     figure;
%     plot(diff(diff(sum_smooth_cwt)))
%     hold on
%     plot((impacttimes-start_time)*Fs_pcb,0,'rx','MarkerSize',8)
%     plot(estimated_impacts(:,1),0,'bx')

end

%% cont from prev section
figure;
plot((impacttimes-start_time)*Fs_pcb,0,'rx','MarkerSize',8)
hold on
plot(estimated_impacts(:,1),0.5,'bx')
ylim([-0.5 1])

% combine impacts
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
    
final_pred
scores
figure;
plot(final_pred,0.5,'bx')
hold on
plot((impacttimes-start_time)*Fs_pcb,0,'rx','MarkerSize',8)
ylim([-0.5 1])






%%
freq_lower = 150; % Hz; frequency limits for cwt
freq_higher = 250;
w=8;
start_time = impacts(segments(w-1),1)/Fs_fsr-0.15;
stop_time = impacts(segments(w)-1,1)/Fs_fsr+0.15;

start_idx_pcb = findTindex(start_time,pcbTime);
stop_idx_pcb = findTindex(stop_time,pcbTime);
    
pkprom = [0.0015,0.0015,0.0015,0.0015,0.0015,0.0015];
pkheight = [0.006,0.016,0.0014,0.004,0.0127,0.009];
for s = 1:sensorN % for each sensor
    pcbclip = wien_pcbD(start_idx_pcb:stop_idx_pcb,s);
    [wt,f] = cwt(pcbclip,Fs_pcb); % uses default Morse wavelet
    valid_f_idx = find(f < freq_higher & f > freq_lower);
    cwt_freq = f(valid_f_idx);
    cwt_mag = abs(wt(valid_f_idx,:));
    sum_cwt = sum(cwt_mag,1);
    sum_smooth_cwt = movmean(sum_cwt, 800);
%     [pks, locs, ~,~] = findpeaks(sum_smooth_cwt, 'MinPeakProminence',0.0015);
    [pks, locs, ~,~] = findpeaks(sum_smooth_cwt,'MinPeakProminence',pkprom(s),'MinPeakHeight',pkheight(s));
    figure;
%     findpeaks(sum_smooth_cwt, 'MinPeakProminence',pkprom(s),'MinPeakHeight',pkheight(s))
    findpeaks(sum_smooth_cwt,'MinPeakProminence',pkprom(s),'MinPeakHeight',pkheight(s))
    hold on
    impacttimes = impacts(segments(w-1):segments(w)-1,1)/Fs_fsr;
    plot((impacttimes-start_time)*Fs_pcb,0,'rx','MarkerSize',8)
end

%% interpolation for overlapping impacts (without wiener) 5/17/22

overlap_idx = find(est_overlapping == 1);
for i = 5:length(overlap_idx) % 3 bc first segment is trash
    curr_idx = overlap_idx(i);
    curr_label = correct_estimates(curr_idx,2);
    if curr_label == 1 | curr_label == 2
        curr_label = 1;
    else
        curr_label = 2;
    end
    % find last impact with same label
    prev_idx = curr_idx - 1;
    prev_label = correct_estimates(prev_idx,2);
    if prev_label == 1 | prev_label == 2
        prev_label = 1;
    else
        prev_label = 2;
    end
    while prev_label ~= curr_label | est_overlapping(prev_idx) == 1
        prev_idx = prev_idx - 1;
        if prev_idx > 0
            prev_label = correct_estimates(prev_idx,2);
        else
            disp("Overlapping at start of impact")
            curr_idx
            break;
        end
    end
    % find next impact with same label
    next_label = correct_estimates(curr_idx+1,2);
    if next_label == 1 | next_label == 2
        next_label = 1;
    else
        next_label = 2;
    end
    next_idx = curr_idx + 1;
    while next_label ~= curr_label | est_overlapping(next_idx) == 1
        next_idx = next_idx + 1;
        if next_idx < correctN
            next_label = correct_estimates(next_idx,2);
        else
            disp("Overlapping at end of impact")
            curr_idx
            break;
        end
    end
    for s = 1:sensorN
        % get current overlapping clip, using same technique as above
        arrival_idx = round(correct_estimates(curr_idx,1)*Fs_pcb);
        est_arrival_idx(curr_idx,s) = arrival_idx;
        next_arrival_idx = curr_idx+1;
        next_arrival = correct_estimates(next_arrival_idx,1)*Fs_pcb;
        while est_arrival_idx(curr_idx,s) == next_arrival
            next_arrival_idx = next_arrival_idx+1;
            next_arrival = correct_estimates(next_arrival_idx,1)*Fs_pcb;
        end
        [lasti,nextimpact] = min([round(correct_estimates(next_arrival_idx,1)*Fs_pcb),arrival_idx + round(0.46*Fs_pcb)]);
        window = wien_pcbD(arrival_idx:lasti,s);
        [wt,f] = cwt(window,Fs_pcb);
        valid_f_idx = find(freq_higher & f > freq_lower);
        cwt_mag = abs(wt(valid_f_idx,:));
        sum_cwt = sum(cwt_mag,1);
        sum_smooth_cwt = movmean(sum_cwt, 800);
        deriv_window = movmean(diff(sum_smooth_cwt),1200);
        [~,minidx] = min(deriv_window);
        [~,zeroidx] = min(abs(deriv_window(minidx:end)));
        abs_last_idx = zeroidx + minidx - 1 + arrival_idx;
        est_last_idx(curr_idx,s) = abs_last_idx;
        
        % extract features w/ interpolate
        est_peak_mag(curr_idx,s) = (est_peak_mag(prev_idx,s) + est_peak_mag(next_idx,s))/2;
        est_energy(curr_idx,s) = (est_energy(prev_idx,s) + est_energy(next_idx,s))/2;
        est_cwt_peak(curr_idx,s) = (est_cwt_peak(prev_idx,s) + est_cwt_peak(next_idx,s))/2;
        est_cwt_energy(curr_idx,s) = (est_cwt_energy(prev_idx,s) + est_cwt_energy(next_idx,s))/2;
    end
    % extract ground truth coordinates
    curr_time = correct_estimates(curr_idx,1);
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
    est_coordinates(curr_idx,:) = [allmocap(mocap_idx,1,mocap_label), allmocap(mocap_idx,3,mocap_label),mocap_label];
end



%%
Fs = 12800;
[impactN, sensorN] = size(arrival_idx);
% looks at the next 3 impacts to label them
real_ID = zeros(impactN,1);
person1_idx = find(impacts(:,4) == 1 | impacts(:,4) == 2);
real_ID(person1_idx) = 1; % label person 1 as 1
person2_idx = find(impacts(:,4) == 3 | impacts(:,4) == 4);
real_ID(person2_idx) = 2; % label person 2 as 2

ID_estimates = zeros(impactN,sensorN);
ID_estimates(1,:) = 2; % always label 1st impact as person 1
ID_estimates(2,:) = 1; % always label 2nd impact as person 2
ID_estimates(3,:) = 2;
ID_confidence = zeros(impactN,sensorN);

for s = 1:sensorN
    for i = 1:impactN
        curr_ID = ID_estimates(i,s);
        curr_clip = filt_pcbD(findTindex(impact_starts(i),pcbTime):findTindex(impact_ends(i,s),pcbTime),s);
        first_clip = [];
        sec_clip = [];
        third_clip = [];
        commonF = linspace(0,2000,1000);
        
        if i+1 <= impactN && ID_estimates(i+1,s) == 0
            first_clip = filt_pcbD(findTindex(impact_starts(i+1),pcbTime):findTindex(impact_ends(i+1,s),pcbTime),s);
        end
        if i+2 <= impactN && ID_estimates(i+2,s) == 0
            sec_clip = filt_pcbD(findTindex(impact_starts(i+2),pcbTime):findTindex(impact_ends(i+2,s),pcbTime),s);
        end
        if i+3 <= impactN && ID_estimates(i+3,s) == 0
            third_clip = filt_pcbD(findTindex(impact_starts(i+3),pcbTime):findTindex(impact_ends(i+3,s),pcbTime),s);
        end
        
        if ~isempty([first_clip;sec_clip;third_clip])
            max_amplitude = [max(abs(curr_clip)),max(abs(first_clip)),max(abs(sec_clip)),max(abs(third_clip))];
            amp_ratios = [];
            for k = 2:length(max_amplitude)
                if ~isempty(max_amplitude(k))
                    amp_ratios(end+1) = max_amplitude(k)/max_amplitude(1);
                else
                    amp_ratios(end+1) = NaN;
                end
            end
            if i < 15
                amp_ratios
            end
            [val,idx] = min(abs(amp_ratios - 1)); % get the idx with ratio closest to 1
            ID_estimates(i+idx,s) = curr_ID;
            % confidence calculates how different it is from max value
            ID_confidence(i+idx,s) = abs(val - max(abs(amp_ratios - 1)));
        end
    end
end

final_estimates = zeros(impactN,1);
% make sensors vote for the most likely ID
% take vote of the lowest uncertainty
for i = 1:impactN
    curr_estimates = ID_estimates(i,:);
    curr_confidence = ID_confidence(i,:);
%     final_estimates(i) = mode(curr_estimates); % vote most often 
    [~,max_id] = max(curr_confidence); % using confidence
    final_estimates(i) = curr_estimates(max_id);
end

figure;
plot(person1_idx,zeros(length(person1_idx),1),'rx')
hold on
plot(person2_idx,zeros(length(person2_idx),1),'bo')
est_person1_idx = find(final_estimates == 1);
est_person2_idx = find(final_estimates == 2);
plot(est_person1_idx,ones(length(est_person1_idx),1),'rx')
plot(est_person2_idx,ones(length(est_person2_idx),1),'bo')
legend('Real person 1','Real person 2','Estimated person 1','Estimated person 2')
ylim([-1 2])

figure;
results = real_ID - final_estimates;
plot(results,'.')
title('Real ID - estimated ID')
ylim([-2 2])
correct_idx = find(results == 0);
accuracy = length(correct_idx)/length(results)

N = length(arrival_idx(:,1));
down_arrival_idx = round(arrival_idx(:,1)./10);
for i = 1:10
    % just sensor 1
    clip = filt_pcbD(arrival_idx(i):arrival_idx(i)+2560,1); % 0.5s window
    clip_d = iddata(clip,[],0.000078125); % sample time is 0.5
    na = 20;
    nc = 20;
%     sys = armax(clip_d,[na nc]);
    sys = ar(clip_d, na);
    figure; compare(clip_d, sys)
end

%% for experiments 12/13/21-12/16/21 initial file to run

clear all
close all
filepath = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\';
subj = '1'; 
% ################# change ########################
intervention = 'regular1'; % regular1 brace1 brace2 weight1 weight2 regular2
MCvalueSet = [3 4 5 6 7 8]; % order of Mocap left/right based on keyset
FSRvalueSet = [1 2]; % order of FSR inputs based on keyset
FSRkeySet = {'Lheel','Rheel'};
% #################################################

filepath = [filepath, 'Subj ',subj,'\subj',subj','_'];
% load 3 datasets
load([filepath, intervention])
load([filepath, 'fsr_', intervention])
T = readtable([filepath, 'mocap_', intervention]);

% constants
num_sensors = 4;
Fs_pcb = 12800;
Mfsr = containers.Map(FSRkeySet, FSRvalueSet);
MCkeySet = {'Lx','Ly','Lz','Rx','Ry','Rz'};
Mmocap = containers.Map(MCkeySet, MCvalueSet);

[mocapT, mocapL, mocapR] = convertMocap(T, Mmocap);

% for i = 239144:266129 % for error in fsr trigger
%     Data(9,i) = -3;
% end

%% data clipping to mocap length 8/30/21
[accData, accTime, fsrData, fsrTime] = clip_accelfsr_fromMocap(Data, Time, Fs);
[pcbData, pcbTime] = clip_pcb_fromMocap(datas, times);

% check this is correct. These numbers should be nearly the same:
mocapT(end)
fsrTime(end)
pcbTime(end)


%% extracting data 9/9/21
% overall plot for visual check
% plot_3data(pcbData,pcbTime,fsrData,fsrTime,mocapR,mocapL,mocapT,Mfsr)

% clean and filter pcb data
downData = downsample(pcbData, 10);
filt_pcbD = lpf_data(pcbData);
filt_down_pcbD = lpf_data(downData);

% finding footfalls based on fsr heel data
% [impacts, Rheel,- Lheel, Rtoe, Ltoe] = findimpacts_fsr(fsrTime,fsrData,Mfsr);
L_dist = 180; % min distance between heel peaks
R_dist = 180;
min_threshL = 25; % value between peaks threshold, if not lower than this then omit peak
min_threshR = 25;

% sometimes the fsr values are 0
for i = 1:length(fsrData)
    if fsrData(i,1) < 10
        fsrData(i,1) = fsrData(i-1,1);
    end
    if fsrData(i,2) < 1
        fsrData(i,2) = fsrData(i-1,2);
    end
end

impacts = findimpacts_fsr_accel(fsrTime,fsrData,Mfsr,L_dist,R_dist,min_threshL,min_threshR);

%% fix small errors in impacts 11/5/21
heel_start_wrong = [27739,37126]; % these need to be same length
heel_start_right = [27688,36986];

heel_pk_wrong = []; % index, these need to be same length
heel_pk_right = [];

impacts = manual_fix_fsr(impacts,fsrData,Mfsr,heel_start_wrong,heel_start_right,heel_pk_wrong,heel_pk_right);

%% to fix subj4_regular1 bc fsr data was bad
% old_impacts = impacts;
% keepimpacts = [];
% for i = 1:length(impacts)
%     if impacts(i,1) < 14085
%         keepimpacts(end+1) = i;
%     elseif impacts(i,1) > 24344
%         keepimpacts(end+1) = i;
%     end
% end
% impacts = old_impacts(keepimpacts,:);
% 
% figure;
% subplot(2,1,1)
% hold on
% plot(fsrData(:,Mfsr('Lheel')))
% lefts = find(impacts(:,4) == 0);
% plot(impacts(lefts,1),fsrData(impacts(lefts,1),Mfsr('Lheel')),'rx','MarkerSize',12)
% plot(impacts(lefts,2),fsrData(impacts(lefts,2),Mfsr('Lheel')),'bx','MarkerSize',12)
% title('Left heel')
% 
% subplot(2,1,2)
% hold on
% plot(fsrData(:,Mfsr('Rheel')))
% rights = find(impacts(:,4) == 1);
% plot(impacts(rights,1),fsrData(impacts(rights,1),Mfsr('Rheel')),'rx','MarkerSize',12)
% plot(impacts(rights,2),fsrData(impacts(rights,2),Mfsr('Rheel')),'bx','MarkerSize',12)
% title('Right heel')

%% extract turning points from looking at x cord 12/13/21

% eliminate any impacts off the floor
x_min = -3500;
x_max = 3600;
elim_i = [1]; % first step is always on the resonance

for i = 1:length(impacts)
    impact_time = fsrTime(impacts(i,1));
    mocap_i = findTindex(impact_time,mocapT);
    if impacts(i,4) == 1 % right foot
        coord = mocapR(mocap_i);
        if coord > x_max || coord < x_min
            elim_i(end+1) = i;
        end
    else % left foot
        coord = mocapL(mocap_i);
        if coord > x_max || coord < x_min
            elim_i(end+1) = i;
        end
    end
end

impacts(elim_i,:) = [];

% visually check
% fsr
RorL = length(impacts(1,:)); % always the last column
figure;
subplot(2,1,1)
hold on
plot(fsrData(:,Mfsr('Lheel')))
lefts = find(impacts(:,RorL) == 0);
plot(impacts(lefts,1),fsrData(impacts(lefts,1),Mfsr('Lheel')),'rx','MarkerSize',12)
plot(impacts(lefts,2),fsrData(impacts(lefts,2),Mfsr('Lheel')),'bx','MarkerSize',12)
title('Left heel')

subplot(2,1,2)
hold on
plot(fsrData(:,Mfsr('Rheel')))
rights = find(impacts(:,RorL) == 1);
plot(impacts(rights,1),fsrData(impacts(rights,1),Mfsr('Rheel')),'rx','MarkerSize',8)
plot(impacts(rights,2),fsrData(impacts(rights,2),Mfsr('Rheel')),'bx','MarkerSize',12)
title('Right heel')

% pcb
figure;
plot(pcbTime,filt_pcbD(:,1))
hold on
plot(fsrTime(impacts(:,1)),0,'r.','MarkerSize',12)

%% manual eliminate turns bc mocap doesn't work

% % eliminate any impacts off the floor
% starts_elim = [3.08578,10.646,17.9223,25.1478,33.1033,40.7816,48.5384,56.3171,63.9938,71.0226,79.1123,86.9313,95.4166,103.019,111.287,119.404,127.321,135.015,142.831,151.43];
% ends_elim = [7.43836,15.1642,22.3268,29.7812,37.9017,45.4147,53.2272,61.1672,68.7294,75.8789,84.0383,91.6504,99.5757,107.9,116.212,124.351,132.146,139.729,147.674,156.196];
% 
% fsr_idx = findTindex(starts_elim(1),fsrTime);
% elim_idx = find(impacts(:,1) < fsr_idx);
% impacts(elim_idx,:) = [];
% 
% for i = 2:length(starts_elim)
%     idx_start = findTindex(starts_elim(i),fsrTime);
%     idx_end = findTindex(ends_elim(i-1),fsrTime);
%     
%     elim_idx = find(impacts(:,1) > idx_end & impacts(:,1) < idx_start);
%     impacts(elim_idx,:) = [];
% end
% 
% fsr_idx = findTindex(ends_elim(end),fsrTime);
% elim_idx = find(impacts(:,1) > fsr_idx);
% impacts(elim_idx,:) = [];
% 
% % visually check
% % fsr
% RorL = length(impacts(1,:)); % always the last column
% figure;
% subplot(2,1,1)
% hold on
% plot(fsrData(:,Mfsr('Lheel')))
% lefts = find(impacts(:,RorL) == 0);
% plot(impacts(lefts,1),fsrData(impacts(lefts,1),Mfsr('Lheel')),'rx','MarkerSize',12)
% plot(impacts(lefts,2),fsrData(impacts(lefts,2),Mfsr('Lheel')),'bx','MarkerSize',12)
% title('Left heel')
% 
% subplot(2,1,2)
% hold on
% plot(fsrData(:,Mfsr('Rheel')))
% rights = find(impacts(:,RorL) == 1);
% plot(impacts(rights,1),fsrData(impacts(rights,1),Mfsr('Rheel')),'rx','MarkerSize',12)
% plot(impacts(rights,2),fsrData(impacts(rights,2),Mfsr('Rheel')),'bx','MarkerSize',12)
% title('Right heel')
% 
% % pcb
% figure;
% plot(pcbTime,filt_pcbD(:,1))
% hold on
% plot(fsrTime(impacts(:,1)),0,'r.','MarkerSize',12)
%% pcb adjust params to capture whole impact 11/5/21
% look at where the impacts are and adjust window_width and offset
% parameters in the next section

close all
Fs = 12800;
secs = 20; % plot 10 secs of data
figure;
for i = 1:4
    subplot(4,1,i)
    plot(pcbTime(1:secs*Fs), filt_pcbD(1:secs*Fs,i))
    hold on
    indeces = find(fsrTime(impacts(:,1)) < secs);
    plot(fsrTime(impacts(indeces,1)),0,'r.','MarkerSize',10)
end

%% pcb and mocap extract

window_width = Fs*0.3; % shaking lasts < N seconds
offset = 0.1; % start window N behind fsr start
[arrival_idx, peak_idx, peak_mag] = findimpacts_pcb(window_width,offset,impacts,fsrTime,pcbTime,filt_pcbD,num_sensors,true);

% mocap extract
[extracted_pts_R, extracted_pts_L, coordinates, whichfoot] = findimpacts_mocap(impacts,fsrTime,mocapT,mocapR,mocapL,true);

%% fix small errors in pcb 11/6/21
arrival_wrong = [13.7388,13.7409,13.7386,13.7367]; % these are matrixes, 4 cols for sensorN and each row is a different wrong impact
arrival_right = [13.6563,13.6555,13.6635,13.6618];
peak_idx_right = [13.6748,13.668,13.6674,13.697];

[arrival_idx,peak_idx,peak_mag] = manual_fix_pcb(arrival_wrong,arrival_right,peak_idx_right,pcbTime,arrival_idx,peak_idx,peak_mag,filt_pcbD);


%% fix small errors in mocap 11/5/21
r_wrong = [24475]; % these need to be same length
r_right = [24713];

l_wrong = []; % index, these need to be same length
l_right = [];

[extracted_pts_R,extracted_pts_L,coordinates] = manual_fix_mocap(r_wrong,r_right,l_wrong,l_right,mocapR,mocapL,extracted_pts_R,extracted_pts_L,coordinates,whichfoot);

%% extract accel data 11/22/21
% gets the peak of abs value of accelerometer values for x,y,z directions
% 4th col in variable 'acc_pks' indicates whether it is LorR leg.
window = 0.3; % seconds
offset = 0.25;

[acc_pks,acc_pk_idx] = find_accel_impacts(impacts,fsrTime,window,offset,accTime,accData,fsrData,Mfsr);

%% saving data to matlab and excel 9/9/21
processedfilepath = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\';
filename = [processedfilepath, 'subj',subj,'_',intervention];

save(filename,'pcbTime','filt_pcbD','arrival_idx','peak_idx','peak_mag', ...
    'fsrTime','fsrData','impacts', ...
    'acc_pks','acc_pk_idx','accTime','accData',...
    'mocapT','mocapR','mocapL','extracted_pts_R','extracted_pts_L','coordinates','whichfoot')
disp(append("Saved as ", filename))

%%
tic
s=5;
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
[diff_array, label_array, estimateID] = notrecursive_stepID_limp_whole(curr_ID_labels,real_impact_times,step_times,step_times_scores,step_o1,step_o2,step_x1,step_x2)
% [diff_array, label_array, estimateID] = optimized_stepID(curr_ID_labels,real_impact_times,step_times,step_o1,step_o2,step_x1,step_x2)
toc