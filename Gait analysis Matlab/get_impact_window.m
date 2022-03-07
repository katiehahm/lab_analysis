function impact_ends = get_impact_window(impact_starts, filt_pcbD, pcbTime)
% returns sensorN x impactN matrix of impact end times
% uses when signal goes below threshold within a defined window
% 3/6/22

filtsize = size(filt_pcbD);
impact_ends = zeros(length(impact_starts),filtsize(2)); % impactN x sensorN
noise_thresh = zeros(filtsize(2),1); % a noise thresh for each sensor
Fs = 12800;
window_length = 0.4; % 0.4 s for max duration of impact

for i = 1:filtsize(2)
    first_clip = filt_pcbD(1:Fs,i);
    noise_thresh(i) = max(abs(first_clip));
end

for s = 1:filtsize(2)
    for i = 1:length(impact_starts)
        curr_noise = noise_thresh(s);
        curr_time = impact_starts(i);
        curr_idx = findTindex(curr_time, pcbTime);
        if curr_time + window_length > pcbTime(end) % last impact
            next_idx = length(pcbTime);
        else
            next_idx = curr_idx + round(window_length*Fs);
        end
        window = filt_pcbD(curr_idx:next_idx,s);
        [up,~] = envelope(window,300,'peak');
        [~,locs] = findpeaks(up); % start searching for below noise thresh after 1st local maxima
        indeces = find(up(locs(1):end) < curr_noise);
        while isempty(indeces) % keep raising threshold until find the signal end
            curr_noise = curr_noise + 0.00001;
            indeces = find(up(locs(1):end) < curr_noise);
        end
        impact_ends(i,s) = pcbTime(indeces(1) + curr_idx + locs(1));
    end
end

end

