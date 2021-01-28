function [res_1, res_2, res_3, res_4] = combine_data_multiplier_ben(data, sensor_code)

    [sens_1_multiplier, sens_2_multiplier, sens_3_multiplier, sens_4_multiplier] = pos_multiplier_ben(data, sensor_code);
    [sens_1_res, sens_2_res, sens_3_res, sens_4_res] = data_organizer(data, sensor_code);
    
    res_1 = horzcat(sens_1_res, sens_1_multiplier);
    res_2 = horzcat(sens_2_res, sens_2_multiplier);
    res_3 = horzcat(sens_3_res, sens_3_multiplier);
    res_4 = horzcat(sens_4_res, sens_4_multiplier);
    
end

