% to be used after experiment4_footfalldetection.m
%% gets arrival idx, peak mag, energy, coord from estimated step times 4/6/22

est_impactN = length(estimated_impacts);
est_coordinates = zeros(est_impactN,2);
est_arrival_idx = zeros(est_impactN,sensorN);
est_peak_idx = zeros(est_impactN,sensorN);
est_peak_mag = zeros(est_impactN,sensorN);
est_energy = zeros(est_impactN,sensorN);
est_cwt_energy = zeros(est_impactN,sensorN);
window_thresh = 0.0825; % seconds of how long an impact lasts
start_idx_thresh = 1000; % count 1000 idx before curr_time to do aicpick
est_overlapping = zeros(est_impactN,1); % value = 1 if two impacts overlap
freq_lower = 100; % Hz; frequency limits for cwt
freq_higher = 350;
est_walk_edges = zeros(est_impactN,1);

for i = 2:est_impactN
    curr_time = estimated_impacts(i,1);
    past_time = estimated_impacts(i-1,1);
    if curr_time == past_time
        est_overlapping(i) = 1;
        est_overlapping(i-1) = 1;
    end
    
end

person1_idx = find(estimated_impacts(:,2) == 1);
person2_idx = find(estimated_impacts(:,2) == 2);
est_walk_edges(person1_idx(1)) = -1;
est_walk_edges(person2_idx(1)) = -1;

for i = 2:length(person1_idx)
    curr_i = person1_idx(i);
    prev_i = person1_idx(i-1);
    time_diff = estimated_impacts(curr_i,1) - estimated_impacts(prev_i,1);
    if time_diff > 1.5
        est_walk_edges(curr_i) = -1;
        est_walk_edges(prev_i) = 1;
    end
end

for i = 2:length(person2_idx)
    curr_i = person2_idx(i);
    prev_i = person2_idx(i-1);
    time_diff = estimated_impacts(curr_i,1) - estimated_impacts(prev_i,1);
    if time_diff > 1.5
        est_walk_edges(curr_i) = -1;
        est_walk_edges(prev_i) = 1;
    end
end

for i = 1:est_impactN
    if est_overlapping(i) ~= 1 % if not an overlapping impact
        for s = 1:sensorN
            % define window end
            curr_time = estimated_impacts(i,1);
            curr_label = estimated_impacts(i,2);
            if i ~= est_impactN % if it isn't the last impact
                next_time = estimated_impacts(i+1,1);
                window_end = min(next_time, curr_time+window_thresh);
            else
                window_end = min(curr_time + window_thresh, length(pcbData)/Fs_pcb);
            end

            % define window start
            curr_idx = curr_time*Fs_pcb;
            start_idx = round(curr_idx - start_idx_thresh);
            curr_window = filt_pcbD(start_idx:round(window_end*Fs_pcb),s);
            
            est_arrival_idx(i,s) = aic_pick(curr_window, 'to_peak')+start_idx;
            curr_window = filt_pcbD(est_arrival_idx(i,s):round(window_end*Fs_pcb),s);
            
            [maxval, maxidx] = max(curr_window);
            est_peak_mag(i,s) = maxval;
            est_peak_idx(i,s) = maxidx + start_idx - 1;
            
            % energy is from start to window_thresh (this may overlap with
            % some of the next impacts check ############################
            energy_window = filt_pcbD(start_idx:start_idx + (window_thresh*Fs_pcb),s);
            % this code is copied from experiment3_analysis.m:
            est_energy(i,s) = sum(abs(energy_window).^2);
            
            % cwt energy
            [wt,f] = cwt(energy_window, Fs_pcb);
            valid_f_idx = find(freq_higher & f > freq_lower);
            cwt_mag = abs(wt(valid_f_idx,:));
            sum_cwt = sum(cwt_mag,1);
            sum_smooth_cwt = movmean(sum_cwt, 800);
            est_cwt_energy(i,s) = sum(abs(sum_smooth_cwt));
            
        end
        
        % define coordinate from closest real impact time
        
        Limpact_idx = find(impacts(:,4) == ((curr_label-1)*2 + 1));
        Rimpact_idx = find(impacts(:,4) == ((curr_label-1)*2 + 2));
        [Lminval,Lcorresponding_real_idx] = min(abs(impacts(Limpact_idx,1)/Fs_fsr - curr_time));
        [Rminval,Rcorresponding_real_idx] = min(abs(impacts(Rimpact_idx,1)/Fs_fsr - curr_time));
        mocap_idx = round(Fs_mocap*(curr_time + 0.05)); % delay for the heel to get "settled"
        if Lminval < Rminval
            mocap_label = (curr_label - 1)*2 + 2;
        else
            mocap_label = (curr_label - 1)*2 + 1;
        end
        est_coordinates(i,:) = [allmocap(mocap_idx,1,mocap_label), allmocap(mocap_idx,3,mocap_label)];
        
    end
end

%% sanity check
figure;
plot(allmocap(:,3,1))
hold on
plabel = find(estimated_impacts(:,2) == 1);
plot(estimated_impacts(plabel,1)*Fs_mocap,est_coordinates(plabel,2),'rx')

for s = 1:6
    figure;
    plot(pcbTime, filt_pcbD(:,s))
    hold on
    plot(est_arrival_idx(:,s)./Fs_pcb,0,'rx')
    plot(est_peak_idx(:,s)./Fs_pcb,est_peak_mag(:,s),'kx')
end



%% save data

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Praneeth experiment 3_11_22\ProcessedData\both_towards1_estimatedStepTimes';
save(data_root_katie,'est_impactN','est_coordinates','est_arrival_idx',...
    'est_peak_idx','est_peak_mag','est_energy',...
    'est_cwt_energy','est_overlapping','est_walk_edges','-append')


%% uses features to create tracking localization csv for each person 4/7/22

for person = 1:2
    featureV = zeros(1,2+6*7+1);
    
    person_idx = find(estimated_impacts(:,2) == person);

    Nsegments = find(est_walk_edges(person_idx) == -1); % find how many segments there are
    trainN = round(length(Nsegments)*0.2); % # segments in one test set (5 fold)
    segN = 0;

    % for the rest of the impacts
    for i = 1:est_impactN
        if estimated_impacts(i,2) == person
            
            min_arr = min(est_arrival_idx(i,:));
            arrivals = est_arrival_idx(i,:)-min_arr;
            xcoord = est_coordinates(i,2)/1000;
            curr_mag = est_peak_mag(i,:);
            curr_energy = est_energy(i,:);
            curr_cwt_energy = est_cwt_energy(i,:);
            

            if est_overlapping(i) ~= 1 % if not overlapping segment
                if est_walk_edges(i) == -1 % start of segment
                    feature = [0, xcoord, curr_mag, arrivals, curr_energy, curr_cwt_energy, zeros(1,19)];
                    featureV(end+1,:) = feature;
                    segN = segN + 1;
                else
                    if xcoord > -3.65 % filter out data under this bc mocap is wrong
                        if mod(segN,trainN) == 0
                            trainValue = floor(segN/trainN);
                        else
                            trainValue = floor(segN/trainN)+1;
                        end

                        if trainValue > 5
                            trainValue = 5;
                        end

                        prev_x = featureV(end,2);
                        prev_mag = featureV(end,3:8);
                        prev_energy = featureV(end,15:20);
                        prev_cwt_energy = featureV(end,21:26);

                        feature = [trainValue,xcoord,curr_mag,arrivals,curr_energy,curr_cwt_energy,curr_mag./prev_mag,curr_energy./prev_energy,curr_cwt_energy./prev_cwt_energy,prev_x];
                        featureV(end+1,:) = feature;
                    end
                end
            end
        end
    end

    featureV(1,:) = []; % initialization

    filename_root = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Praneeth experiment 3_11_22\ProcessedData\both_towards1_localization_person';
    filename = [filename_root, int2str(person),'.csv'];
    writematrix(featureV,filename)
end
% then run recursive_localization2pp.py

%% try wiener filter to localize overlapping impacts too 4/12/22

% overlapping impacts numbers 13/14
figure; 
overlap = 13;
overlap_idx = estimated_impacts(overlap,1)*Fs_pcb;
plot(pcbTime(overlap_idx - Fs_pcb:overlap_idx + Fs_pcb),filt_pcbD(overlap_idx - Fs_pcb:overlap_idx + Fs_pcb,1))
hold on
plot(estimated_impacts(overlap,1),0,'rx')
plot(estimated_impacts(overlap+1,1),0,'bx')

persono = find(estimated_impacts(overlap-3:overlap-1,2) == 1);
personx = find(estimated_impacts(overlap-3:overlap-1,2) == 2);

persono_idx = persono(end)+overlap-4;
personx_idx = personx(end)+overlap-4;

plot(estimated_impacts(persono_idx,1),0,'ro')
plot(estimated_impacts(personx_idx,1),0,'rx')
window_thresh_idx = round(window_thresh*Fs_pcb);

overlap_clip = filt_pcbD(round(5.15719*Fs_pcb):round(5.15719*Fs_pcb)+window_thresh_idx,1);
persono_clip = filt_pcbD(est_arrival_idx(persono_idx,1):est_arrival_idx(persono_idx,1)+window_thresh_idx,1);
personx_clip = filt_pcbD(est_arrival_idx(personx_idx,1):est_arrival_idx(personx_idx,1)+window_thresh_idx,1);

[yhat_o, H_o] = wienerFilter(persono_clip,overlap_clip,1,true,Fs_pcb);
[yhat_x, H_x] = wienerFilter(personx_clip,overlap_clip,1,true,Fs_pcb);

figure; plot(yhat_o)
figure; plot(yhat_x)

%% same as section 1, but add overlapping segments w/ wiener filter 4/12/22

est_impactN = length(estimated_impacts);
est_coordinates = zeros(est_impactN,2);
est_arrival_idx = zeros(est_impactN,sensorN);
est_peak_idx = zeros(est_impactN,sensorN);
est_peak_mag = zeros(est_impactN,sensorN);
est_energy = zeros(est_impactN,sensorN);
est_cwt_energy = zeros(est_impactN,sensorN);
window_thresh = 0.0825; % seconds of how long an impact lasts
start_idx_thresh = 1000; % count 1000 idx before curr_time to do aicpick
est_overlapping = zeros(est_impactN,1); % value = 1 if two impacts overlap
freq_lower = 100; % Hz; frequency limits for cwt
freq_higher = 350;
est_walk_edges = zeros(est_impactN,1);

for i = 2:est_impactN
    curr_time = estimated_impacts(i,1);
    past_time = estimated_impacts(i-1,1);
    if curr_time == past_time
        est_overlapping(i) = 1;
        est_overlapping(i-1) = 1;
    end
    
end

person1_idx = find(estimated_impacts(:,2) == 1);
person2_idx = find(estimated_impacts(:,2) == 2);
est_walk_edges(person1_idx(1)) = -1;
est_walk_edges(person2_idx(1)) = -1;

for i = 2:length(person1_idx)
    curr_i = person1_idx(i);
    prev_i = person1_idx(i-1);
    time_diff = estimated_impacts(curr_i,1) - estimated_impacts(prev_i,1);
    if time_diff > 1.5
        est_walk_edges(curr_i) = -1;
        est_walk_edges(prev_i) = 1;
    end
end

for i = 2:length(person2_idx)
    curr_i = person2_idx(i);
    prev_i = person2_idx(i-1);
    time_diff = estimated_impacts(curr_i,1) - estimated_impacts(prev_i,1);
    if time_diff > 1.5
        est_walk_edges(curr_i) = -1;
        est_walk_edges(prev_i) = 1;
    end
end

for i = 1:est_impactN
    if est_overlapping(i) ~= 1 % if not an overlapping impact
        for s = 1:sensorN
            % define window end
            curr_time = estimated_impacts(i,1);
            curr_label = estimated_impacts(i,2);
            if i ~= est_impactN % if it isn't the last impact
                next_time = estimated_impacts(i+1,1);
                window_end = min(next_time, curr_time+window_thresh);
            else
                window_end = min(curr_time + window_thresh, length(pcbData)/Fs_pcb);
            end

            % define window start
            curr_idx = curr_time*Fs_pcb;
            start_idx = round(curr_idx - start_idx_thresh);
            curr_window = filt_pcbD(start_idx:round(window_end*Fs_pcb),s);
            
            est_arrival_idx(i,s) = aic_pick(curr_window, 'to_peak')+start_idx;
            curr_window = filt_pcbD(est_arrival_idx(i,s):round(window_end*Fs_pcb),s);
            
            [maxval, maxidx] = max(curr_window);
            est_peak_mag(i,s) = maxval;
            est_peak_idx(i,s) = maxidx + start_idx - 1;
            
            % energy is from start to window_thresh (this may overlap with
            % some of the next impacts check ############################
            energy_window = filt_pcbD(start_idx:start_idx + (window_thresh*Fs_pcb),s);
            % this code is copied from experiment3_analysis.m:
            est_energy(i,s) = sum(abs(energy_window).^2);
            
            % cwt energy
            [wt,f] = cwt(energy_window, Fs_pcb);
            valid_f_idx = find(freq_higher & f > freq_lower);
            cwt_mag = abs(wt(valid_f_idx,:));
            sum_cwt = sum(cwt_mag,1);
            sum_smooth_cwt = movmean(sum_cwt, 800);
            est_cwt_energy(i,s) = sum(abs(sum_smooth_cwt));
            
        end
        
        % define coordinate from closest real impact time
        
        Limpact_idx = find(impacts(:,4) == ((curr_label-1)*2 + 1));
        Rimpact_idx = find(impacts(:,4) == ((curr_label-1)*2 + 2));
        [Lminval,Lcorresponding_real_idx] = min(abs(impacts(Limpact_idx,1)/Fs_fsr - curr_time));
        [Rminval,Rcorresponding_real_idx] = min(abs(impacts(Rimpact_idx,1)/Fs_fsr - curr_time));
        mocap_idx = round(Fs_mocap*(curr_time + 0.05)); % delay for the heel to get "settled"
        if Lminval < Rminval
            mocap_label = (curr_label - 1)*2 + 2;
        else
            mocap_label = (curr_label - 1)*2 + 1;
        end
        est_coordinates(i,:) = [allmocap(mocap_idx,1,mocap_label), allmocap(mocap_idx,3,mocap_label)];
    
    else % overlapping segment
        
        
    end
end

for refernece
    figure; 
overlap = 13;
overlap_idx = estimated_impacts(overlap,1)*Fs_pcb;
plot(pcbTime(overlap_idx - Fs_pcb:overlap_idx + Fs_pcb),filt_pcbD(overlap_idx - Fs_pcb:overlap_idx + Fs_pcb,1))
hold on
plot(estimated_impacts(overlap,1),0,'rx')
plot(estimated_impacts(overlap+1,1),0,'bx')

persono = find(estimated_impacts(overlap-3:overlap-1,2) == 1);
personx = find(estimated_impacts(overlap-3:overlap-1,2) == 2);

persono_idx = persono(end)+overlap-4;
personx_idx = personx(end)+overlap-4;

plot(estimated_impacts(persono_idx,1),0,'ro')
plot(estimated_impacts(personx_idx,1),0,'rx')
window_thresh_idx = round(window_thresh*Fs_pcb);

overlap_clip = filt_pcbD(round(5.15719*Fs_pcb):round(5.15719*Fs_pcb)+window_thresh_idx,1);
persono_clip = filt_pcbD(est_arrival_idx(persono_idx,1):est_arrival_idx(persono_idx,1)+window_thresh_idx,1);
personx_clip = filt_pcbD(est_arrival_idx(personx_idx,1):est_arrival_idx(personx_idx,1)+window_thresh_idx,1);

[yhat_o, H_o] = wienerFilter(persono_clip,overlap_clip,1,true,Fs_pcb);
[yhat_x, H_x] = wienerFilter(personx_clip,overlap_clip,1,true,Fs_pcb);

figure; plot(yhat_o)
figure; plot(yhat_x)






































