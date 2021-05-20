function [joined_steps, dst] = ProcessMoCap(motion_data)
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
    
    stop_idx = 0;

    for i=1:length(x)
        if isnan(x(i))
            x(i) = 0;
        end
        if isnan(y(i))
            y(i) = 0;
        end
        if isnan(z(i))
            z(i) = 0;
        end
        if z(i) > 350
            z(i) = 350;
        end
        if z(i) < 30 && z(i) ~= 0
            if stop_idx == 0 
                stop_idx = i;
            end
        end
        if isnan(t(i))
            t(i) = 0;
        end
    end

    [~, idx_l] = findpeaks(-lz(1:stop_idx),'MinPeakDistance', 200);
    [~, idx_r] = findpeaks(-rz(1:stop_idx),'MinPeakDistance', 200);
    step_l = zeros(length(idx_l),3);
    step_r = zeros(length(idx_r),3);
    
    for s = 1:length(idx_l)
        step_l_idx = idx_l(s);
        step_l(s, :) = [lx(step_l_idx) -ly(step_l_idx) step_l_idx];
    end
    for l = 1:length(idx_r)
        step_r_idx = idx_r(l);
        step_r(l, :) = [rx(step_r_idx) -ry(step_r_idx) step_r_idx];
    end

    sensor_x = [-3585 -3606 3616 3607];
    sensor_y = [3395 -196 -212 3376];

    joined_steps = sortrows(vertcat(step_r, step_l),3);

    dst_1 = zeros(length(joined_steps),1);
    dst_2 = zeros(length(joined_steps),1);
    dst_3 = zeros(length(joined_steps),1);
    dst_4 = zeros(length(joined_steps),1);

    for r = 1:length(joined_steps)
        dst_1(r) = power((power((sensor_x(1) - joined_steps(r, 1)),2) + power((sensor_y(1) - joined_steps(r, 2)),2)),1/2);
        dst_2(r) = power((power((sensor_x(2) - joined_steps(r, 1)),2) + power((sensor_y(2) - joined_steps(r, 2)),2)),1/2);
        dst_3(r) = power((power((sensor_x(3) - joined_steps(r, 1)),2) + power((sensor_y(3) - joined_steps(r, 2)),2)),1/2);
        dst_4(r) = power((power((sensor_x(4) - joined_steps(r, 1)),2) + power((sensor_y(4) - joined_steps(r, 2)),2)),1/2);
    end

    dst = horzcat(dst_1, dst_2, dst_3, dst_4);
    
    figure();
    hold on
    scatter(lx,-ly,10, lz);
    scatter(rx,-ry,10, rz);
    scatter(step_l(:,1), step_l(:, 2), 75, 'MarkerEdgeColor', 'magenta',... 
                                            'MarkerFaceColor', 'none',...
                                            'LineWidth', 2);
    scatter(step_r(:,1), step_r(:, 2), 75, 'MarkerEdgeColor', 'red',... 
                                            'MarkerFaceColor', 'none',...
                                            'LineWidth', 2);
    scatter(sensor_x, sensor_y, 75, 'MarkerEdgeColor', 'black',... 
                                    'MarkerFaceColor', 'black',...
                                    'LineWidth', 2);
    xlabel('x');
    ylabel('y');
    c = colorbar;
    c.Label.String = 'z';
    grid on;
end