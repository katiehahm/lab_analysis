function [onset_idx, peak_idx, peak_val] = TDOA_edited(clean_data,impactN)
    d = clean_data;
    sn = length(d(1,:)); % number of sensors
    d_length = length(d(:,1)); % number of datapoints
    onset_idx = zeros(sn,impactN);
    peak_idx = zeros(sn,impactN);
    peak_val = zeros(sn,impactN);

    % ADJUSTABLE THRESHOLD VALUES
    onset2pk = 500; % max num of datapoints between onset and peak
    onset_thresh = 0.002; % Volts. threshold where noise becomes impact
    impactWidth = 7500; % impact lasts N datapoints

    % find max recursively
    for i = 1:sn
        for j = 1:impactN
            [maxval, maxidx] = max(d(:,i));
            peak_val(i,j) = maxval;
            peak_idx(i,j) = maxidx;
            min_idx = maxidx - onset2pk;
            max_idx = min_idx + impactWidth;
            d(min_idx:max_idx,i) = 0;
        end
    end

    % sorting first to last impact
    [peak_idx,bidx] = sort(peak_idx,2);
    for i = 1:sn
        peak_val(i,:) = peak_val(i,bidx(i,:));
    end
    for i = 1:sn
        prev = peak_idx(i,1);
        for j = 2:length(peak_idx(1,:))
            curr = peak_idx(i,j);
            if (curr - prev) < impactWidth/2 % found 2 peaks in 1 impact
                [maxval, maxidx] = max(d(:,i));
                peak_val(i,j) = maxval;
                peak_idx(i,j) = maxidx;
                min_idx = maxidx - onset2pk;
                max_idx = min_idx + impactWidth;
                d(min_idx:max_idx,i) = 0;
            end
            prev = curr;
        end
    end

    % for each peak, find corresponding onset index
    d = clean_data;
    for i = 1:sn
        for j = 1:impactN
            pki = peak_idx(i,j);
            chunk2 = d(pki-onset2pk:pki,i);
            v = find(chunk2 > onset_thresh);
            if isempty(v)
                [~,oidx] = max(chunk2);
                onset_idx(i,j) = oidx + pki - onset2pk;
            else
                onset_idx(i,j) = v(1) + pki-onset2pk;
            end
        end
    end

    if any(onset_idx==0,'all')
        disp("Not all onsets found")
    end
    if any(peak_idx==0,'all')
        disp("Not all peaks found")
    end
end