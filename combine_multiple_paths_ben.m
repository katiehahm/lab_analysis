function [sensor_1, sensor_2, sensor_3, sensor_4] = combine_multiple_paths_ben(datas, codes)
    sensor_1 = NaN;
    sensor_2 = NaN;
    sensor_3 = NaN;
    sensor_4 = NaN;
    
    for i = 1:length(datas)
        data = datas{i};
        code = codes(i);
        
        if isnan(sensor_1)
            [sensor_1, sensor_2, sensor_3, sensor_4] = combine_data_multiplier_ben(data, code);
        else
            [update_1, update_2, update_3, update_4] = combine_data_multiplier_ben(data, code);
            sensor_1 = vertcat(sensor_1, update_1);
            sensor_2 = vertcat(sensor_2, update_2);
            sensor_3 = vertcat(sensor_3, update_3);
            sensor_4 = vertcat(sensor_4, update_4);            
        end  
    end
    
    figure
    scatter(sensor_1(:,2), sensor_1(:,1))
    improvePlot_ben
    grid on
    title("Sensor 1")
    ylabel("Impact Magnitude")
    xlabel("Distance (Inches)")
    
    figure
    scatter(sensor_2(:,2), sensor_2(:,1))
    improvePlot_ben
    grid on
    title("Sensor 2")
    ylabel("Impact Magnitude")
    xlabel("Distance (Inches)")
    
    figure
    scatter(sensor_3(:,2), sensor_3(:,1))
    improvePlot_ben
    grid on
    title("Sensor 3")
    ylabel("Impact Magnitude")
    xlabel("Distance (Inches)")
    
    figure
    scatter(sensor_4(:,2), sensor_4(:,1))
    improvePlot_ben
    grid on
    title("Sensor 4")
    ylabel("Impact Magnitude")
    xlabel("Distance (Inches)")
    
end
