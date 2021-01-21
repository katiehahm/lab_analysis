function [sens_1_res, sens_2_res, sens_3_res, sens_4_res, sens_length] = data_organizer(data, sensor_code)
    sens_1 = data(:, 1);
    sens_2 = data(:, 2);
    sens_3 = data(:, 3);
    sens_4 = data(:, 4);
    
    [peaks_1, idx_1, width_1, prominence_1] = TDOA3(sens_1);
    [peaks_2, idx_2, width_2, prominence_2] = TDOA3(sens_2);
    [peaks_3, idx_3, width_3, prominence_3] = TDOA3(sens_3);
    [peaks_4, idx_4, width_4, prominence_4] = TDOA3(sens_4);
    
    if sensor_code == 11100000000
        pos_1_peak_1 = peaks_1(1:5,1);
        pos_1_peak_2 = peaks_2(1:5,1);
        pos_1_peak_3 = peaks_3(1:5,1);
        pos_1_peak_4 = peaks_4(1:5,1);
        
        pos_2_peak_1 = peaks_1(6:10,1);
        pos_2_peak_2 = peaks_2(6:10,1);
        pos_2_peak_3 = peaks_3(6:10,1);
        pos_2_peak_4 = peaks_4(6:10,1);
        
        pos_3_peak_1 = peaks_1(11:15,1);
        pos_3_peak_2 = peaks_2(11:15,1);
        pos_3_peak_3 = peaks_3(11:15,1);
        pos_3_peak_4 = peaks_4(11:15,1);
        
        sens_1_res = vertcat(pos_1_peak_1, pos_2_peak_1, pos_3_peak_1);
        sens_2_res = vertcat(pos_1_peak_2, pos_2_peak_2, pos_3_peak_2);
        sens_3_res = vertcat(pos_1_peak_3, pos_2_peak_3, pos_3_peak_3);
        sens_4_res = vertcat(pos_1_peak_4, pos_2_peak_4, pos_3_peak_4);
    end
    if sensor_code == 00011100000
        pos_1_peak_1 = peaks_1(1:5,1);
        pos_1_peak_2 = peaks_2(1:5,1);
        pos_1_peak_3 = peaks_3(1:5,1);
        pos_1_peak_4 = peaks_4(1:5,1);
        
        pos_2_peak_1 = peaks_1(6:10,1);
        pos_2_peak_2 = peaks_2(6:10,1);
        pos_2_peak_3 = peaks_3(6:10,1);
        pos_2_peak_4 = peaks_4(6:10,1);
        
        pos_3_peak_1 = peaks_1(11:15,1);
        pos_3_peak_2 = peaks_2(11:15,1);
        pos_3_peak_3 = peaks_3(11:15,1);
        pos_3_peak_4 = peaks_4(11:15,1);
        
        sens_1_res = vertcat(pos_1_peak_1, pos_2_peak_1, pos_3_peak_1);
        sens_2_res = vertcat(pos_1_peak_2, pos_2_peak_2, pos_3_peak_2);
        sens_3_res = vertcat(pos_1_peak_3, pos_2_peak_3, pos_3_peak_3);
        sens_4_res = vertcat(pos_1_peak_4, pos_2_peak_4, pos_3_peak_4);
    end
    if sensor_code == 00000011100
        pos_1_peak_1 = peaks_1(1:5,1);
        pos_1_peak_2 = peaks_2(1:5,1);
        pos_1_peak_3 = peaks_3(1:5,1);
        pos_1_peak_4 = peaks_4(1:5,1);
        
        pos_2_peak_1 = peaks_1(6:10,1);
        pos_2_peak_2 = peaks_2(6:10,1);
        pos_2_peak_3 = peaks_3(6:10,1);
        pos_2_peak_4 = peaks_4(6:10,1);
        
        pos_3_peak_1 = peaks_1(11:15,1);
        pos_3_peak_2 = peaks_2(11:15,1);
        pos_3_peak_3 = peaks_3(11:15,1);
        pos_3_peak_4 = peaks_4(11:15,1);
        
        sens_1_res = vertcat(pos_1_peak_1, pos_2_peak_1, pos_3_peak_1);
        sens_2_res = vertcat(pos_1_peak_2, pos_2_peak_2, pos_3_peak_2);
        sens_3_res = vertcat(pos_1_peak_3, pos_2_peak_3, pos_3_peak_3);
        sens_4_res = vertcat(pos_1_peak_4, pos_2_peak_4, pos_3_peak_4);
    end
    if sensor_code == 00000000011
        pos_1_peak_1 = peaks_1(1:5,1);
        pos_1_peak_2 = peaks_2(1:5,1);
        pos_1_peak_3 = peaks_3(1:5,1);
        pos_1_peak_4 = peaks_4(1:5,1);
        
        pos_2_peak_1 = peaks_1(6:10,1);
        pos_2_peak_2 = peaks_2(6:10,1);
        pos_2_peak_3 = peaks_3(6:10,1);
        pos_2_peak_4 = peaks_4(6:10,1);
        
        sens_1_res = vertcat(pos_1_peak_1, pos_2_peak_1);
        sens_2_res = vertcat(pos_1_peak_2, pos_2_peak_2);
        sens_3_res = vertcat(pos_1_peak_3, pos_2_peak_3);
        sens_4_res = vertcat(pos_1_peak_4, pos_2_peak_4);
    end
    sens_length = length(sens_1_res);
end