% Simplified code for TDOA
% This allows us to take in data without cleaning it.
% We also do not need to include the number of impacts.
% Feel free to adjust thresholds, just comment out originals.
function [peaks, idx, width, prominence] = TDOA3(data)
%     filtered_data = lpf_data(data.datas)
    min_height = 0.015;                   % default min height
    min_distance = 5000;                % default min distance between peaks
    min_prominence = 0.015;                 % default min peak prominence
    [peaks, idx, width, prominence] = findpeaks(data(:,1),... 
                                                'MinPeakHeight', min_height,... 
                                                'MinPeakDistance', min_distance,... 
                                                'MinPeakProminence', min_prominence);
end 