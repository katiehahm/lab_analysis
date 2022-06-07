% for use after experiment4_allprocessingcompiled.m
% gets features and makes csv for localization estimation
%% feature extraction w/ cwt approach 5/14/22

clear diff % clear this as a variable name
% constants
est_impactN = length(clean_final_estimates_labels);
window_thresh = 0.0825; % seconds of how long an impact lasts
start_idx_thresh = 1000; % count 1000 idx before curr_time to do aicpick
freq_lower = 100; % Hz; frequency limits for cwt
freq_higher = 350;
% noise_thresh = [0.000285,0.00026,0.00018,0.00026,0.00024,0.00056]; % determined empirically by looking at plot

correct_estimates = [0,0,0];
correct_coords = [0,0,0];
correct_ta = [0];
% find all correct impacts & ground truth coords
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
correct_estimates(1,:) = [];
correct_coords(1,:) = [];
correct_ta(1,:) = [];
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

segments_list = unique(correct_estimates(:,3));
segmentsN = segments_list(end)-segments_list(1) + 1;
% extract values
for i = 1:length(segments_list)
    seg = segments_list(i);
    seg_idx = find(correct_estimates(:,3) == (seg));
    seg_estimates = correct_estimates(seg_idx,:);
    estimatesN = length(seg_estimates(:,1));
    for s = 1:sensorN
        % extract cwt of segment
        not_overlap = [];
        first_idx = round((seg_estimates(1,1)-0.5)*Fs_pcb);
        last_idx = min([round((seg_estimates(estimatesN,1)+0.5)*Fs_pcb),length(wien_pcbD(:,1))]);
        seg_pcb = wien_pcbD(first_idx:last_idx,s);
        [wt,f] = cwt(seg_pcb,Fs_pcb);
        valid_f_idx = find(freq_higher & f > freq_lower);
        cwt_mag = abs(wt(valid_f_idx,:));
        sum_cwt = sum(cwt_mag,1);
        sum_smooth_cwt = movmean(sum_cwt, 800);
        pcb2cwt_idx = first_idx; % add this to result to get pcb idx
        % extract features from cwt sum
        for i = 1:estimatesN
            % check for overlapping
            if i~= estimatesN & round(seg_estimates(i,1),3) == round(seg_estimates(i+1,1),3)
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
                    window2 = wien_pcbD(arrival_idx:abs_last_idx,s);
                    [wt2,f2] = cwt(window2,Fs_pcb);
                    valid_f_idx2 = find(freq_higher & f2 > freq_lower);
                    cwt_mag2 = abs(wt2(valid_f_idx2,:));
                    sum_cwt2 = sum(cwt_mag2,1);
                    sum_smooth_cwt2 = movmean(sum_cwt2, 800);
                    est_cwt_energy(seg_idx(i),s) = sum(sum_smooth_cwt2.^2);
                    % peak mag
                    est_peak_mag(seg_idx(i),s) = max(abs(window2));
                    % peak cwt mag
                    est_cwt_peak(seg_idx(i),s) = max(sum_smooth_cwt2);
                    % energy
                    est_energy(seg_idx(i),s) = sum(abs(window2).^2);
                end
                % ground truth coordinates
                curr_time = correct_estimates(seg_idx(i),1);
                curr_label = correct_estimates(seg_idx(i),2);
                if curr_label == 1 | curr_label == 2
                    curr_label = 1;
                else
                    curr_label = 2;
                end
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

%% interpolation for overlapping impacts (without wiener) 5/17/22

overlap_idx = find(est_overlapping == 1);
for i = 1:length(overlap_idx) 
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
        next_arrival = round(correct_estimates(next_arrival_idx,1)*Fs_pcb);
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


%% sanity check 5/13/22

% to check if overlapping detection is correct
% red x's should be vertically stacked
figure;
plot(correct_estimates(:,1),linspace(1,correctN,correctN),'b.')
hold on
plot(correct_estimates(find(est_overlapping == 1),1),find(est_overlapping == 1),'rx')

% to test that missing features are only due to overlapping impacts
% should run bc the arrays are same length
% arrival_zeros = find(est_arrival_idx(:,1) == 0);
% figure;
% plot(arrival_zeros,find(est_overlapping == 1))

% to check if impact window is accurate
arrival_nonzeros = find(est_arrival_idx(:,1) ~= 0);
for s = 1:sensorN
    figure;
    plot(pcbTime, wien_pcbD(:,s))
    hold on
    plot(est_arrival_idx(:,s)./Fs_pcb,0,'rx','MarkerSize',14)
    plot(est_last_idx(:,s)./Fs_pcb,0,'gx','MarkerSize',9)
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

% check if peak & energy extraction are correct
% zoom in to the signal, ignore the zero and negative values
% peaks should match
test_start = 70;
test_end = 95;
test_seg = wien_pcbD(est_arrival_idx(test_start,1):est_arrival_idx(test_end,1),1);
[wt,f] = cwt(test_seg,Fs_pcb);
valid_f_idx = find(freq_higher & f > freq_lower);
cwt_mag = abs(wt(valid_f_idx,:));
sum_cwt = sum(cwt_mag,1);
sum_smooth_cwt = movmean(sum_cwt, 800);
figure; plot(sum_smooth_cwt)
hold on
plot(est_arrival_idx(test_start:test_end,1)-est_arrival_idx(test_start,1),est_cwt_peak(test_start:test_end,1),'rx')
est_cwt_energy(test_start:test_end,1)
figure;
plot(test_seg)
hold on
plot(est_arrival_idx(test_start:test_end,1)-est_arrival_idx(test_start,1),est_peak_mag(test_start:test_end,1),'rx')
est_energy(test_start:test_end,1)

%% save features to mat file
save(filename,'correct_coords','correct_estimates','est_arrival_idx','est_cwt_energy','est_cwt_peak','est_energy',...
    'est_overlapping','est_peak_mag','est_last_idx','correct_ta','correctN','-append')
disp(append("Saved as ", filename))

%% make localization csv with overlapping impacts for tracking 5/16/22

featureV_p1 = zeros(1,53);
featureV_p2 = zeros(1,53);

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

% move files from ProcessedData to ExcelData
% then run exp4_recursive_localization.py OR
% exp4_notrecursive_localization.py

% %% make TA csv with overlapping impacts 5/18/22
% 
% % copy over excel files, delete first row and column!
% filename = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Jenny 1\ProcessedData\ExcelData\both_regular2_localization_p2_results.csv';
% T = readtable(filename);
% A = table2array(T);
% real_loc = A(:,1);
% predict_loc = A(:,2);
% 
% filename = 'C:/Users/Katie/Dropbox (MIT)/Lab/Analysis/Experiment4/Jenny 1/ProcessedData/ExcelData/both_regular2_localization_p2_withta.csv';
% T = readtable(filename);
% A = table2array(T);
% ta_arr = A(:,4);
% features = A(:,5:28);
% 
% s1 = [-3.590,-3.343];
% s2 = [-3.580,2.61];
% s3 = [3.639,2.11];
% s4 = [3.650,-3.412];
% % these are values from interpolating s1-4
% s5 = [3.61,2.36];
% s6 = [3.62,-3.3775];
% featureV = zeros(1,31);
% 
% for i = 1:length(real_loc)
%     xcoord = predict_loc(i);
%     % only using x coord
%     % ############################################################## check!
%     dist1 = sqrt( (xcoord-s1(1)).^2 );
%     dist2 = sqrt( (xcoord-s2(1)).^2 );
%     dist3 = sqrt( (xcoord-s3(1)).^2 );
%     dist4 = sqrt( (xcoord-s4(1)).^2 );
%     dist5 = sqrt( (xcoord-s5(1)).^2 );
%     dist6 = sqrt( (xcoord-s6(1)).^2 );
% %     dist1 = sqrt( (xcoord-s1(1)).^2 + (ycoord-s1(2)).^2 );
% %     dist2 = sqrt( (xcoord-s2(1)).^2 + (ycoord-s2(2)).^2 );
% %     dist3 = sqrt( (xcoord-s3(1)).^2 + (ycoord-s3(2)).^2 );
% %     dist4 = sqrt( (xcoord-s4(1)).^2 + (ycoord-s4(2)).^2 );
% %     dist5 = sqrt( (xcoord-s5(1)).^2 + (ycoord-s5(2)).^2 );
% %     dist6 = sqrt( (xcoord-s6(1)).^2 + (ycoord-s6(2)).^2 );
% 
%     % acc, dist 1-6, cwt peak, cwt energy, peak, energy
%     feature = [abs(ta_arr(i)),dist1,dist2,dist3,dist4,dist5,dist6,features(i,:)];
%     featureV(end+1,:) = feature;
% end
% featureV(1,:) = [];
% filename = 'C:/Users/Katie/Dropbox (MIT)/Lab/Analysis/Experiment4/Jenny 1/ProcessedData/ExcelData/both_regular2_ta_p2.csv';
% writematrix(featureV,filename)
% 
% % then run exp4_TAestimation.py

%% making a localization csv for all takes 5/24/22

takes = {'regular1', 'limp1', 'limp2', 'weight1', 'weight2', 'regular2'};
newA = zeros(1,53);
person = '2';
for t = 1:length(takes)
    filename = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\ExcelData\both_', char(takes(t)),'_localization_p',person','_withta.csv'];
    T = readtable(filename);
    A = table2array(T);
    newA = [newA;A];
end
newA(1,:) = []; % initialization

filename = ['C:/Users/Katie/Dropbox (MIT)/Lab/Analysis/Experiment4/April 3/ProcessedData/ExcelData/alltakes_localization_p',person','_withta.csv'];
writematrix(newA,filename)

%% make TA csv with overlapping impacts 5/26/22 using alltakes localization results

filename1 = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\ExcelData\alltakes_localization_p1_results.csv';
T = readtable(filename1);
loc_results1 = table2array(T);
filename2 = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\ExcelData\alltakes_localization_p2_results.csv';
T = readtable(filename2);
loc_results2 = table2array(T);

takes = {'regular1', 'limp1', 'limp2', 'weight1', 'weight2', 'regular2'};
s1 = [-3.590,-3.343];
s2 = [-3.580,2.61];
s3 = [3.639,2.11];
s4 = [3.650,-3.412];
% these are values from interpolating s1-4
s5 = [3.61,2.36];
s6 = [3.62,-3.3775];


for p = 1:2
    person = int2str(p);
    rowcount = 1;
    for t = 1:length(takes)
        filename = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\ExcelData\both_', char(takes(t)),'_localization_p',person','_withta.csv'];
        T = readtable(filename);
        A = table2array(T);
        numrows = length(A(:,1));
        if p == 1
            pred_locs = loc_results1(:,2); % second col has est values
        else
            pred_locs = loc_results2(:,2);
        end
        locs = pred_locs(rowcount:rowcount + numrows - 1);
        featureV = zeros(1,31);
        for i = 1:numrows
            xcoord = locs(i);
            dist1 = sqrt( (xcoord-s1(1)).^2 );
            dist2 = sqrt( (xcoord-s2(1)).^2 );
            dist3 = sqrt( (xcoord-s3(1)).^2 );
            dist4 = sqrt( (xcoord-s4(1)).^2 );
            dist5 = sqrt( (xcoord-s5(1)).^2 );
            dist6 = sqrt( (xcoord-s6(1)).^2 );
            feature = [A(i,4),dist1,dist2,dist3,dist4,dist5,dist6,A(i,5:28)];
            featureV(end+1,:) = feature;
        end
        featureV(1,:) = []; % initialization
        newfilename = ['C:/Users/Katie/Dropbox (MIT)/Lab/Analysis/Experiment4/April 3/ProcessedData/ExcelData/both_', char(takes(t)),'_ta_p',person,'.csv'];
        writematrix(featureV,newfilename)
        rowcount = rowcount + numrows;
    end
end

% then run exp4_TAestimation.py

%% doing kmeans on TA estimation results all takes 5/24/22

takes = {'regular1', 'limp1', 'limp2', 'weight1', 'weight2', 'regular2'};
newA = zeros(1,53);
TA_rmse = 0;
figure;
hold on
for t = 1:length(takes)
    filename1 = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\ExcelData\both_', char(takes(t)),'_ta_p1_results.csv'];
    T1 = readtable(filename1);
    A1 = table2array(T1);
    est_TA1 = A1(2:end,3);
    [~,cent1] = kmeans(est_TA1,2);
    
    filename2 = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\ExcelData\both_', char(takes(t)),'_ta_p2_results.csv'];
    T2 = readtable(filename2);
    A2 = table2array(T2);
    est_TA2 = A2(2:end,3);
    [~,cent2] = kmeans(est_TA2,2);
    
    load(['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\both_',char(takes(t))])
    foot_labels = correct_coords(:,3);
    p1f1_idx = find(foot_labels == 1);
    p1f2_idx = find(foot_labels == 2);
    p2f1_idx = find(foot_labels == 3);
    p2f2_idx = find(foot_labels == 4);
    real_p1f1 = mean(correct_ta(p1f1_idx));
    real_p1f2 = mean(correct_ta(p1f2_idx));
    real_p2f1 = mean(correct_ta(p2f1_idx));
    real_p2f2 = mean(correct_ta(p2f2_idx));
    
    plot(min(cent1),min([real_p1f1,real_p1f2]),'ro')
    plot(max(cent1),max([real_p1f1,real_p1f2]),'bo')
    plot(min(cent2),min([real_p2f1,real_p2f2]),'rx')
    plot(max(cent2),max([real_p2f1,real_p2f2]),'bx')
    legend('Leg 1 person 1','Leg 2 person 1','Leg 1 person 2','Leg 2 person 2')
    
    % calculate rmse
    TA_rmse = TA_rmse + abs(min(cent1)-min([real_p1f1,real_p1f2]));
    TA_rmse = TA_rmse + abs(max(cent1)-max([real_p1f1,real_p1f2]));
    TA_rmse = TA_rmse + abs(min(cent2)-min([real_p2f1,real_p2f2]));
    TA_rmse = TA_rmse + abs(max(cent2)-max([real_p2f1,real_p2f2]));
end

xlabel('Estimated TA values (g)')
ylabel('Measured TA values (g)')
title('TA estimation performance for each leg across all interventions')
xlim([1 5])
ylim([1 5])
TA_rmse/(4*length(takes)) % final rmse

%% make localization csv without overlapping impacts 5/16/22
% 
% featureV_p1 = zeros(1,26);
% featureV_p2 = zeros(1,26);
% for i = 1:length(est_coordinates)
%     if est_overlapping(i) == 0 % not overlapping
%         feature = [correct_coords(i,1),abs(correct_ta(i)),est_cwt_peak(i,:),est_cwt_energy(i,:),est_peak_mag(i,:),est_energy(i,:)];
%         if correct_estimates(i,2) == 1
%             featureV_p1(end+1,:) = feature;
%         else
%             featureV_p2(end+1,:) = feature;
%         end
%     end
% end
% featureV_p1(1,:) = [];
% featureV_p2(1,:) = [];
% filename1 = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Jenny 1\ProcessedData\ExcelData\both_regular1_localization_p1_withta.csv';
% filename2 = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Jenny 1\ProcessedData\ExcelData\both_regular1_localization_p2_withta.csv';
% writematrix(featureV_p1,filename1)
% writematrix(featureV_p2,filename2)
% 
% %% make TA csv without overlapping impacts 5/16/22
% 
% filename = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Jenny 1\ProcessedData\ExcelData\both_regular1_localization_p1';
% T = readtable([filename,'_results.csv']);
% A = table2array(T);
% real_loc = A(:,1);
% predict_loc = A(:,2);
% 
% T = readtable([filename,'_withta.csv']);
% A = table2array(T);
% ta_arr = A(:,2);
% features = A(:,3:end);
% 
% s1 = [-3.590,-3.343];
% s2 = [-3.580,2.61];
% s3 = [3.639,2.11];
% s4 = [3.650,-3.412];
% % these are values from interpolating s1-4
% s5 = [3.61,2.36];
% s6 = [3.62,-3.3775];
% featureV = zeros(1,31);
% 
% for i = 1:length(real_loc)
%     xcoord = predict_loc(i);
%     % only using x coord
%     % ############################################################## check!
%     dist1 = sqrt( (xcoord-s1(1)).^2 );
%     dist2 = sqrt( (xcoord-s2(1)).^2 );
%     dist3 = sqrt( (xcoord-s3(1)).^2 );
%     dist4 = sqrt( (xcoord-s4(1)).^2 );
%     dist5 = sqrt( (xcoord-s5(1)).^2 );
%     dist6 = sqrt( (xcoord-s6(1)).^2 );
% %     dist1 = sqrt( (xcoord-s1(1)).^2 + (ycoord-s1(2)).^2 );
% %     dist2 = sqrt( (xcoord-s2(1)).^2 + (ycoord-s2(2)).^2 );
% %     dist3 = sqrt( (xcoord-s3(1)).^2 + (ycoord-s3(2)).^2 );
% %     dist4 = sqrt( (xcoord-s4(1)).^2 + (ycoord-s4(2)).^2 );
% %     dist5 = sqrt( (xcoord-s5(1)).^2 + (ycoord-s5(2)).^2 );
% %     dist6 = sqrt( (xcoord-s6(1)).^2 + (ycoord-s6(2)).^2 );
% 
%     % acc, dist 1-6, cwt peak, cwt energy, peak, energy
%     feature = [abs(ta_arr(i)),dist1,dist2,dist3,dist4,dist5,dist6,features(i,:)];
%     featureV(end+1,:) = feature;
% end
% featureV(1,:) = [];
% filename = [filename, '_TA.csv'];
% writematrix(featureV,filename)


%% wiener filter & interpolation for overlapping impacts 5/17/22

% overlap_idx = find(est_overlapping == 1);
% for i = 7:length(overlap_idx) % put 7 here bc first segment is trash
%     curr_idx = overlap_idx(i); % idx within length of clean_final...
%     data_idx = curr_idx*Fs_pcb; % idx within data segment
%     curr_time = correct_estimates(curr_idx,1);
%     curr_label = correct_estimates(curr_idx,2);
%     % find last impact with same label
%     prev_idx = curr_idx - 1;
%     prev_label = correct_estimates(prev_idx,2);
%     while prev_label ~= curr_label | est_overlapping(prev_idx) == 1
%         prev_idx = prev_idx - 1;
%         if prev_idx > 0
%             prev_label = correct_estimates(prev_idx,2);
%         else
%             disp("Overlapping at start of impact")
%             quit()
%         end
%     end
%     % find next impact with same label
%     next_label = correct_estimates(curr_idx+1,2);
%     next_idx = curr_idx + 1;
%     while next_label ~= curr_label | est_overlapping(next_label) == 1
%         next_idx = next_idx + 1;
%         if next_idx < correctN
%             next_label = correct_estimates(next_idx,2);
%         else
%             disp("Overlapping at end of impact")
%             quit()
%         end
%     end
%     for s = 1:sensorN
%         % previous clip
%         prev_clip = filt_pcbD(est_arrival_idx(prev_idx,s):est_last_idx(prev_idx,s),s);
%         % next clip
%         next_clip = filt_pcbD(est_arrival_idx(next_idx,s):est_last_idx(next_idx,s),s);
%         % get current overlapping clip, using same technique as above
%         arrival_idx = round(correct_estimates(curr_idx,1)*Fs_pcb);
%         est_arrival_idx(curr_idx,s) = arrival_idx;
%         next_arrival_idx = curr_idx+1;
%         next_arrival = correct_estimates(next_arrival_idx,1)*Fs_pcb;
%         while est_arrival_idx(curr_idx,s) == next_arrival
%             next_arrival_idx = next_arrival_idx+1;
%             next_arrival = correct_estimates(next_arrival_idx,1)*Fs_pcb;
%         end
%         [lasti,nextimpact] = min([round(correct_estimates(next_arrival_idx,1)*Fs_pcb),arrival_idx + round(0.46*Fs_pcb)]);
%         window = filt_pcbD(arrival_idx:lasti,s);
%         [wt,f] = cwt(window,Fs_pcb);
%         valid_f_idx = find(freq_higher & f > freq_lower);
%         cwt_mag = abs(wt(valid_f_idx,:));
%         sum_cwt = sum(cwt_mag,1);
%         sum_smooth_cwt = movmean(sum_cwt, 800);
%         deriv_window = movmean(diff(sum_smooth_cwt),1200);
%         [~,minidx] = min(deriv_window);
%         [~,zeroidx] = min(abs(deriv_window(minidx:end)));
%         abs_last_idx = zeroidx + minidx - 1 + arrival_idx;
%         est_last_idx(curr_idx,s) = abs_last_idx;
%         overlap_clip = filt_pcbD(arrival_idx:abs_last_idx,s);
%         % make clips same length
%         prev_clipL = length(prev_clip);
%         next_clipL = length(next_clip);
%         overlap_clipL = length(overlap_clip);
%         [~,longclip] = max([prev_clipL,next_clipL,overlap_clipL]);
%         if longclip == 1 % longest is prev clip
%             next_clip = wiener_clipL_adjust(prev_clipL, next_clipL, s, filt_pcbD, est_arrival_idx(next_idx,s), est_last_idx(next_idx,s), noise_thresh(s));
%             overlap_clip = wiener_clipL_adjust(prev_clipL, overlap_clipL, s, filt_pcbD, arrival_idx, abs_last_idx, noise_thresh(s));
%         elseif longclip == 2 % longest is next clip
%             prev_clip = wiener_clipL_adjust(next_clipL, prev_clipL, s, filt_pcbD, est_arrival_idx(prev_idx,s), est_last_idx(prev_idx,s), noise_thresh(s));
%             overlap_clip = wiener_clipL_adjust(next_clipL, overlap_clipL, s, filt_pcbD, arrival_idx, abs_last_idx, noise_thresh(s));
%         else % longest is overlapping clip
%             prev_clip = wiener_clipL_adjust(overlap_clipL, prev_clipL, s, filt_pcbD, est_arrival_idx(prev_idx,s), est_last_idx(prev_idx,s), noise_thresh(s));
%             next_clip = wiener_clipL_adjust(overlap_clipL, next_clipL, s, filt_pcbD, est_arrival_idx(next_idx,s), est_last_idx(next_idx,s), noise_thresh(s));
%         end
%         % perform wiener filter on closer clip
%         if (next_idx - curr_idx) > (curr_idx - prev_idx)
%             [yhat, ~] = wienerFilter(prev_clip,overlap_clip,1,true,Fs_pcb);
%         else
%             [yhat, ~] = wienerFilter(next_clip,overlap_clip,1,true,Fs_pcb);
%         end
%         % interpolate to find scaling factor
%         new_mag = (max(abs(prev_clip)) + max(abs(next_clip)))/2;
%         scaling_factor = new_mag/max(abs(yhat));
%         yhat = yhat.*scaling_factor;
%         % extract energy from extracted signal
%         est_energy(curr_idx,s) = sum(abs(yhat).^2);
%         % extract peak from extracted signal
%         est_peak_mag(curr_idx,s) = new_mag;
%         % extract cwt energy & peak from extracted signal
%         [wt,f] = cwt(yhat, Fs_pcb);
%         valid_f_idx = find(freq_higher & f > freq_lower);
%         cwt_mag = abs(wt(valid_f_idx,:));
%         sum_cwt = sum(cwt_mag,1);
%         sum_smooth_cwt = movmean(sum_cwt, 800);
%         est_cwt_energy(curr_idx,s) = sum(abs(sum_smooth_cwt));
%         est_cwt_peak(curr_idx,s) = max(sum_smooth_cwt);
%     end
% end

































