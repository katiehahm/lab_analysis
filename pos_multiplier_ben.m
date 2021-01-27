function [sens_1_multiplier, sens_2_multiplier, sens_3_multiplier, sens_4_multiplier] = pos_multiplier_ben(data, sensor_code)
    [sens_1_res, sens_2_res, sens_3_res, sens_4_res, len] = data_organizer(data, sensor_code);
    
    sens_1_multiplier = NaN;
    sens_2_multiplier = NaN;
    sens_3_multiplier = NaN;
    sens_4_multiplier = NaN;
    
    if sensor_code == 11100000000
        positions = [1 2 3];
    end
    if sensor_code == 00011100000
        positions = [4 5 6];
    end
    if sensor_code == 00000011100
        positions = [7 8 9];
    end
    if sensor_code == 00000000011
        positions = [10 11];
    end
    
    for pos_num = positions
        if pos_num == 1
            dist_to_1 = 125;
            dist_to_2 = 65;
            dist_to_3 = 77;
            dist_to_4 = 155;
        end
        if pos_num == 2
            dist_to_1 = 115;
            dist_to_2 = 77;
            dist_to_3 = 87;
            dist_to_4 = 122;
        end
        if pos_num == 3
            dist_to_1 = 98;
            dist_to_2 = 91;
            dist_to_3 = 100;
            dist_to_4 = 106;
        end
        if pos_num == 4
            dist_to_1 = 83;
            dist_to_2 = 107;
            dist_to_3 = 115;
            dist_to_4 = 92;
        end
        if pos_num == 5
            dist_to_1 = 70;
            dist_to_2 = 124;
            dist_to_3 = 131;
            dist_to_4 = 81;
        end
        if pos_num == 6
            dist_to_1 = 61;
            dist_to_2 = 142;
            dist_to_3 = 148;
            dist_to_4 = 73;
        end
        if pos_num == 7
            dist_to_1 = 38;
            dist_to_2 = 134;
            dist_to_3 = 162;
            dist_to_4 = 98;
        end
        if pos_num == 8
            dist_to_1 = 51;
            dist_to_2 = 115;
            dist_to_3 = 146;
            dist_to_4 = 104;
        end
        if pos_num == 9
            dist_to_1 = 68;
            dist_to_2 = 96;
            dist_to_3 = 132;
            dist_to_4 = 113;
        end
        if pos_num == 10
            dist_to_1 = 85;
            dist_to_2 = 154;
            dist_to_3 = 138;
            dist_to_4 = 50;60
        end
        if pos_num == 11
            dist_to_1 = 110;
            dist_to_2 = 169;
            dist_to_3 = 132;
            dist_to_4 = 28;
        end
        num_trials = len/length(positions);
        if isnan(sens_1_multiplier)
            sens_1_multiplier = dist_to_1.* ones(num_trials,1);
            sens_2_multiplier = dist_to_2.* ones(num_trials,1);
            sens_3_multiplier = dist_to_3.* ones(num_trials,1);
            sens_4_multiplier = dist_to_4.* ones(num_trials,1);
        else
            sens_1_multiplier = vertcat(sens_1_multiplier, dist_to_1.*ones(num_trials,1));
            sens_2_multiplier = vertcat(sens_2_multiplier, dist_to_2.*ones(num_trials,1));
            sens_3_multiplier = vertcat(sens_3_multiplier, dist_to_3.*ones(num_trials,1));
            sens_4_multiplier = vertcat(sens_4_multiplier, dist_to_4.*ones(num_trials,1));
        end
    end
end

