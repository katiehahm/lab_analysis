function [final] = Prep_MoCap(motion_data, sensor_data)
    [motion_data_crop,sensor_data_crop, ~, ~] = TimeShift(motion_data,sensor_data);

    [joined_steps, dst] = ProcessMoCap(motion_data_crop);
    [sensor_data_crop_clean] = lpf_data(sensor_data_crop);
    [~, peak_idx, peak_val] = TDOA_edited(sensor_data_crop_clean,length(joined_steps));

    sen_1 = horzcat(peak_idx(1,:).', peak_val(1,:).');
    sen_2 = horzcat(peak_idx(2,:).', peak_val(2,:).');
    sen_3 = horzcat(peak_idx(3,:).', peak_val(3,:).');
    sen_4 = horzcat(peak_idx(4,:).', peak_val(4,:).');

    sen_1_sorted = sortrows(sen_1, 1);
    sen_2_sorted = sortrows(sen_2, 1);
    sen_3_sorted = sortrows(sen_3, 1);
    sen_4_sorted = sortrows(sen_4, 1);

    sen_merged_1 = horzcat(sen_1_sorted(:, 1), sen_2_sorted(:, 1), sen_3_sorted(:, 1), sen_4_sorted(:, 1));
    sen_merged_2 = horzcat(sen_1_sorted(:, 2), sen_2_sorted(:, 2), sen_3_sorted(:, 2), sen_4_sorted(:, 2));
    sen_merged = horzcat(sen_merged_1, sen_merged_2);

    final = zeros(length(joined_steps), 12);
    for r = 1:length(joined_steps)
        first_idx = 1000000000000;
        high_mag = 0;
        for c = 1:4
            if sen_merged(r, c) < first_idx
                first_idx = sen_merged(r, c);
            end
        end
        for c = 5:8
            if sen_merged(r, c) > high_mag
                high_mag = sen_merged(r, c);
            end
        end
        for c = 1:4
            final(r, c) = sen_merged(r, c) - first_idx;
        end
        for c = 5:8
            final(r, c) = sen_merged(r, c) - high_mag;
        end
        for c = 9:12
            final(r,c) = dst(r,c-8);
        end
    end
end