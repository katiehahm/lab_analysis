function [motion_data_crop,sensor_data_crop, sensor_start, motion_start] = TimeShift(motion_data,sensor_data)
    start_mag = 0;
    start_idx = 0;
    
    lx = table2array(motion_data(:,9));
    ly = table2array(motion_data(:,11));
    lz = table2array(motion_data(:,10));
    lt = table2array(motion_data(:,2));

    rx = table2array(motion_data(:,18));
    ry = table2array(motion_data(:,20));
    rz = table2array(motion_data(:,19));
    rt = table2array(motion_data(:,2));

    x = vertcat(lx, rx);
    y = vertcat(ly, ry);
    z = vertcat(lz, rz);
    t = vertcat(lt, rt);

    for i=1:length(z)
        if isnan(z(i))
            z(i) = 0;
        end
        if z(i) > start_mag
            start_mag = z(i);
            start_idx = i - length(z)/2;
        end
    end

    [peaks_l, idx_l] = findpeaks(-lz,'MinPeakDistance', 200);
    [peaks_r, idx_r] = findpeaks(-rz,'MinPeakDistance', 200);
    step_l = zeros(length(idx_l),3);
    step_r = zeros(length(idx_r),3);
    
    for s = 1:length(idx_l)
        step_l_idx = idx_l(s);
        step_l(s, :) = [lx(step_l_idx) -ly(step_l_idx) step_l_idx];
    end
    for s = 1:length(idx_r)
        step_r_idx = idx_r(s);
        step_r(s, :) = [rx(step_l_idx) -ry(step_l_idx) step_r_idx];
    end

    motion_start = 0;
    for num = length(step_r):-1:1
        if step_r(num,3) > start_idx
            motion_start = step_r(num, 3);
        end
    end
    
    max_val = 0;
    sensor_start = 0;
    for i = 1:length(sensor_data)
        for c = 1:4
            if sensor_data(i,c) > max_val
                max_val = sensor_data(i,c);
                sensor_start = i;
            end
        end
    end
        
    motion_data_crop = motion_data(motion_start:height(motion_data),:);
    sensor_data_crop = sensor_data(sensor_start-5000:length(sensor_data),:);
end

