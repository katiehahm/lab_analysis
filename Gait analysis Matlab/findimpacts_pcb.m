function [arrival_idx, peak_idx, peak_mag] = findimpacts_pcb(impacts,fsrTime,pcbT,filt_pcbD,Fs,num_sensors, visualize)
% finds the accelerometer impacts based on impacts found through fsr
% uses the fsr timestamps to window region and extract
% arrival index, peak index, and peak mag
% if no impact was found, value is NaN
% 6/8/21


window_width = Fs*0.3; % shaking lasts < 0.3s
offset = 0.25; % start window 1/4 of window behind the fsr start
impactT = fsrTime(impacts(:,1));
arrival_idx = zeros(length(impactT),num_sensors);
peak_idx = zeros(length(impactT),num_sensors);
peak_mag = zeros(length(impactT),num_sensors);
for i = 1:length(impactT)
    pcbi = findTindex(impactT(i),pcbT);
    starti = pcbi - offset*window_width;
    for j = 1:num_sensors
        endi = min(starti+window_width, length(filt_pcbD(:,1)));
        window = filt_pcbD(starti:endi,j);
        if any(window > 0)
            arrival_idx(i,j) = aic_pick(window, 'to_peak')+starti;
            [mag,idx] = max(abs(window));
            peak_mag(i,j) = mag;
            peak_idx(i,j) = idx + starti;
        else
            arrival_idx(i,j) = NaN;
            peak_idx(i,j) = NaN;
            peak_mag(i,j) = NaN;
        end
    end
end

% delete NaN values:
% eliminate = find(all(isnan(arrival_idx),2));
% arrival_idx(eliminate,:) = [];
% peak_idx(eliminate,:) = [];
% peak_mag(eliminate,:) = [];
% impactT(eliminate) = [];
% impacts(eliminate,:) = [];

% visualize
if visualize
    pk_idx = peak_idx;
    pk_mag = peak_mag;
    ar_idx = arrival_idx;
    pk_idx(isnan(peak_idx)) = 1;
    pk_mag(isnan(peak_mag)) = 1;
    ar_idx(isnan(arrival_idx)) = 1;
    figure;
    for i=1:num_sensors
        subplot(4,1,i)
        plot(pcbT,filt_pcbD(:,i))
        hold on
        plot(pcbT(pk_idx(:,i)),pk_mag(:,i),'r.','MarkerSize',10)
        plot(pcbT(ar_idx(:,i)),0,'k.','MarkerSize',10)
    end
end

end

