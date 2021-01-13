function [onset_idx, peak_idx, peak_val] = TDOA_copy(d,impactN,Fs,loc_names)
    sn = length(d(1,:)); % number of sensors
    d_length = length(d(:,1)); % number of datapoints
    % this only works with 3 sensors******
    onset_idx1 = [];
    peak_idx1 = [];
    peak_val1 = [];
    onset_idx2 = [];
    peak_idx2 = [];
    peak_val2 = [];
    onset_idx3 = [];
    peak_idx3 = [];
    peak_val3 = [];

    % number of datapoints after an onset that should not average near thresh
    impactR = round(Fs/300); % impact lasts > 1/nth of a second
    % threshold where value is an impact, not noise
    thresh = 0.008; % Volts
    % number of datapoints after an onset that should include a peak
    pkR = round(Fs/15); % impact lasts > 1/15th of a second
    impactL = round(Fs/2.5); % whole impact length where multiple onsets cannot exist

    s=1;
    for i = 1:(d_length-max(impactR,pkR)-1)
        if d(i,s) ~= 0
            avg_seg = sum(d(i:i+impactR,s))/impactR;
            if (avg_seg > thresh)
                if (length(onset_idx1) == 0) || (i-onset_idx1(end))> impactL
                    onset_idx1(end+1) = i;
                    [pk, pk_i] = max(d(i:i+pkR,s));
                    peak_idx1(end+1) = pk_i+i;
                    peak_val1(end+1) = pk;
                end
            end

        end
    end
    s=2;
    for i = 1:(d_length-max(impactR,pkR)-1)
        if d(i,s) ~= 0
            avg_seg = sum(d(i:i+impactR,s))/impactR;
            if (avg_seg > thresh)
                if (length(onset_idx2) == 0) || (i-onset_idx2(end))> impactL
                    onset_idx2(end+1) = i;
                    [pk, pk_i] = max(d(i:i+pkR,s));
                    peak_idx2(end+1) = pk_i+i;
                    peak_val2(end+1) = pk;
                end
            end
        end
    end
    s=3;
    for i = 1:(d_length-max(impactR,pkR)-1)
        if d(i,s) ~= 0
            avg_seg = sum(d(i:i+impactR,s))/impactR;
            if (avg_seg > thresh)
                if (length(onset_idx3) == 0) || (i-onset_idx3(end))> impactL
                onset_idx3(end+1) = i;
                [pk, pk_i] = max(d(i:i+pkR,s));
                peak_idx3(end+1) = pk_i+i;
                peak_val3(end+1) = pk;
                end
            end
        end
    end

    % for testing******
%     figure;
%     hold on
%     subplot(3,1,1)
%     hold on
%     plot(d(:,1))
%     plot(onset_idx1,0,'bx')
%     plot(peak_idx1,peak_val1,'ro')
%     subplot(3,1,2)
%     hold on
%     plot(d(:,2))
%     plot(onset_idx2,0,'bx')
%     plot(peak_idx2,peak_val2,'ro')
%     subplot(3,1,3)
%     hold on
%     plot(d(:,3))
%     plot(onset_idx3,0,'bx')
%     plot(peak_idx3,peak_val3,'ro')



    % only select impacts that are felt by all 3 sensors
    % ensures all indeces are same length
    timeDelayT = 300; % +-threshold of number of indeces due to time delay 
    onset_idx = zeros(impactN, sn);
    peak_idx = zeros(impactN, sn);
    peak_val = zeros(impactN, sn);
    count = 1;
    % loop through all elements of first array
    % if a similar impact is found by the other two arrays, store it
    for i = 1:length(onset_idx1)
        val = onset_idx1(i);
        range_min = val - timeDelayT;
        range_max = val + timeDelayT;
        [x2,y2] = find(onset_idx2>range_min & onset_idx2<range_max);
        [x3,y3] = find(onset_idx3>range_min & onset_idx3<range_max);
        if ~isempty(x2) && ~isempty(x3)
            if (x2 && x3) % onset idx found within range for both
                onset_idx(count,1) = onset_idx1(i);
                onset_idx(count,2) = onset_idx2(y2);
                onset_idx(count,3) = onset_idx3(y3);
                peak_idx(count,1) = peak_idx1(i);
                peak_idx(count,2) = peak_idx2(y2);
                peak_idx(count,3) = peak_idx3(y3);
                peak_val(count,1) = peak_val1(i);
                peak_val(count,2) = peak_val2(y2);
                peak_val(count,3) = peak_val3(y3);
                count = count + 1;
            end
        end
    end

    if onset_idx(impactN,1) == 0
        disp("Not all impacts found")
    end

    figure;
    subplot(3,1,1)
    title(append('Filtered data at ', loc_names(1)))
    xlabel('Impact number')
    ylabel('Volts (V)')
    hold on
    plot(d(:,1))
    plot(onset_idx(:,1),0,'bx')
    plot(peak_idx(:,1),peak_val(:,1),'ro')
    subplot(3,1,2)
    title(append('Filtered data at ', loc_names(2)))
    xlabel('Impact number')
    ylabel('Volts (V)')
    hold on
    plot(d(:,2))
    plot(onset_idx(:,2),0,'bx')
    plot(peak_idx(:,2),peak_val(:,2),'ro')
    subplot(3,1,3)
    title(append('Filtered data at ', loc_names(3)))
    xlabel('Impact number')
    ylabel('Volts (V)')
    hold on
    plot(d(:,3))
    plot(onset_idx(:,3),0,'bx')
    plot(peak_idx(:,3),peak_val(:,3),'ro')

end

