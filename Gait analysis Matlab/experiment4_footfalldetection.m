% used to detect footfalls using cwt
% 3/24/22

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Praneeth experiment 3_11_22\ProcessedData\both_towards1';
load(string(data_root_katie))

freq_lower = 100; % Hz; frequency limits for cwt
freq_higher = 350;

estimated_impacts = [0,0]; % [impact time, sensor N], for each segment
final_estimates = [0,0]; % [impact time, segment ID num]
impact_thresh = 0.04; % (s) of how close an impact should be to the last

arrival_window = 3000; % 5000/Fs_pcb width of window to find arrival time
impact_time_thresh = 0.1; % if two impacts are 0.1s within each other, combine

for w = 2:2%length(segments) % for each walking segment
    start_time = impacts(segments(w-1)+1,1)/Fs_fsr-1;
    stop_time = impacts(segments(w),1)/Fs_fsr+1;
    
    start_idx_pcb = findTindex(start_time,pcbTime);
    stop_idx_pcb = findTindex(stop_time,pcbTime);
    
    impacttimes = impacts(segments(w-1)+1:segments(w),1)/Fs_fsr;
    
    for s = 1:sensorN % for each sensor
        pcbclip = pcbData(start_idx_pcb:stop_idx_pcb,s);
        [wt,f] = cwt(pcbclip,Fs_pcb); % uses default Morse wavelet
        valid_f_idx = find(f < freq_higher & f > freq_lower);
        cwt_freq = f(valid_f_idx);
        cwt_mag = abs(wt(valid_f_idx,:));
        sum_cwt = sum(cwt_mag,1);
        sum_smooth_cwt = movmean(sum_cwt, 800);
        [pks, locs, ~,~] = findpeaks(sum_smooth_cwt, 'MinPeakProminence',0.0015);
%         figure;
%         findpeaks(sum_smooth_cwt, 'MinPeakProminence',0.0015)
%         hold on
%         impacttimes = impacts(segments(w-1)+1:segments(w),1)/Fs_fsr;
%         plot((impacttimes-start_time)*Fs_pcb,0,'rx','MarkerSize',8)
        for i = 1:length(locs)
            if i == 1
                starti = max(locs(i) - arrival_window, 1);
            else
                [minval, minidx] = min(sum_smooth_cwt(locs(i-1):locs(i)));
                starti = max(locs(i) - arrival_window, locs(i-1)+minidx - arrival_window/5);
            end
            window = sum_smooth_cwt(starti:locs(i));
            arrive_idx = aic_pick(window, 'to_peak')+starti;
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
end

final_estimates(1,:) = []; % from initialization









