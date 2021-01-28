function [sens_1_res, sens_2_res, sens_3_res, sens_4_res, sens_length] = data_organizer(data, sensor_code)

    if sensor_code == 11100000000 || sensor_code == 00011100000 || sensor_code == 00000011100
        pos_1_peak_1 = data(1:5,1);
        pos_1_peak_2 = data(1:5,2);
        pos_1_peak_3 = data(1:5,3);
        pos_1_peak_4 = data(1:5,4);
        
        pos_2_peak_1 = data(6:10,1);
        pos_2_peak_2 = data(6:10,2);
        pos_2_peak_3 = data(6:10,3);
        pos_2_peak_4 = data(6:10,4);
        
        pos_3_peak_1 = data(11:15,1);
        pos_3_peak_2 = data(11:15,2);
        pos_3_peak_3 = data(11:15,3);
        pos_3_peak_4 = data(11:15,4);
        
        sens_1_res = vertcat(pos_1_peak_1, pos_2_peak_1, pos_3_peak_1);
        sens_2_res = vertcat(pos_1_peak_2, pos_2_peak_2, pos_3_peak_2);
        sens_3_res = vertcat(pos_1_peak_3, pos_2_peak_3, pos_3_peak_3);
        sens_4_res = vertcat(pos_1_peak_4, pos_2_peak_4, pos_3_peak_4);
    end

    if sensor_code == 00000000011
        pos_1_peak_1 = data(1:5,1);
        pos_1_peak_2 = data(1:5,1);
        pos_1_peak_3 = data(1:5,1);
        pos_1_peak_4 = data(1:5,1);
        
        pos_2_peak_1 = data(6:10,1);
        pos_2_peak_2 = data(6:10,1);
        pos_2_peak_3 = data(6:10,1);
        pos_2_peak_4 = data(6:10,1);
        
        sens_1_res = vertcat(pos_1_peak_1, pos_2_peak_1);
        sens_2_res = vertcat(pos_1_peak_2, pos_2_peak_2);
        sens_3_res = vertcat(pos_1_peak_3, pos_2_peak_3);
        sens_4_res = vertcat(pos_1_peak_4, pos_2_peak_4);
    end
    sens_length = length(sens_1_res);
end