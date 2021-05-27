% Created by: Ben (bendwyer@mit.edu) 5/21/21

% Inputs: 
%   motion_data - data loaded from motion capture system 
%   sensor_data - data loaded from sensor system (.datas)

% Purpose: 
%   This function is meant to combine all the other helper functions 
%   into one, easy to use function. We only need to input our data and
%   it will produce the array we need to input into localization. This
%   code will not produce any intermediate values.

% Output:
%   final - array that contains the values needed to run localization.

function [final] = Prep_MoCap(motion_data, sensor_data)
    [motion_data_crop,sensor_data_crop, ~, ~] = TimeShift(motion_data,sensor_data);

    [joined_steps, dst] = ProcessMoCap(motion_data_crop);
    [sensor_data_crop_clean] = lpf_data(sensor_data_crop);
    [onset_idx, peak_idx, peak_val] = TDOA_edited(sensor_data_crop_clean,length(joined_steps));

    sen_1 = horzcat(onset_idx(1,:).', peak_idx(1,:).', peak_val(1,:).');
    sen_2 = horzcat(onset_idx(2,:).', peak_idx(2,:).', peak_val(2,:).');
    sen_3 = horzcat(onset_idx(3,:).', peak_idx(3,:).', peak_val(3,:).');
    sen_4 = horzcat(onset_idx(4,:).', peak_idx(4,:).', peak_val(4,:).');

    sen_1_sorted = sortrows(sen_1, 2);
    sen_2_sorted = sortrows(sen_2, 2);
    sen_3_sorted = sortrows(sen_3, 2);
    sen_4_sorted = sortrows(sen_4, 2);

    sen_merged_1 = horzcat(sen_1_sorted(:, 1), sen_2_sorted(:, 1), sen_3_sorted(:, 1), sen_4_sorted(:, 1));
    sen_merged_2 = horzcat(sen_1_sorted(:, 2), sen_2_sorted(:, 2), sen_3_sorted(:, 2), sen_4_sorted(:, 2));
    sen_merged_3 = horzcat(sen_1_sorted(:, 3), sen_2_sorted(:, 3), sen_3_sorted(:, 3), sen_4_sorted(:, 3));
    sen_merged = horzcat(sen_merged_1, sen_merged_2, sen_merged_3);

    final = zeros(length(joined_steps), 16);
    for r = 1:length(joined_steps)
        first_peak_idx = 1000000000000;
        first_onset_idx = 1000000000000;
        high_mag = 0;
        for c = 1:4
            if sen_merged(r, c) < first_peak_idx
                first_peak_idx = sen_merged(r, c);
            end
        end
        for c = 5:8
            if sen_merged(r, c) < first_onset_idx
                first_onset_idx = sen_merged(r, c);
            end
        end
        for c = 9:12
            if sen_merged(r, c) > high_mag
                high_mag = sen_merged(r, c);
            end
        end
        for c = 1:4
            final(r, c) = sen_merged(r, c) - first_onset_idx;
        end
        for c = 5:8
            final(r, c) = sen_merged(r, c) - first_peak_idx;
        end
        for c = 9:12
            final(r, c) = sen_merged(r, c) - high_mag;
        end
        for c = 13:16
            final(r,c) = dst(r,c-12);
        end
    end
    label = ["onset_idx_diff_1", "onset_idx_diff_2", "onset_idx_diff_3", "onset_idx_diff_4",...
             "peak_idx_diff_1", "peak_idx_diff_2", "peak_idx_diff_3", "peak_idx_diff_4",...
             "mag_diff_1", "mag_diff_2", "mag_diff_3", "mag_diff_4",...
             "dst_to_1", "dst_to_2", "dst_to_3", "dst_to_4", "x", "y"];
    final = horzcat(final, joined_steps(:,1:2));
    final = vertcat(label,final);
end