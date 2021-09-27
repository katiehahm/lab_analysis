function [distance] = distance_btw_coordinates(coord_idx1, coord_idx2,mocap, axis)
% calculates the distance traveled by the foot between two heel strikes
% used to find any outliers (turns) to extract walking segments
% coord_idx1 < coord_idx2 always
% axis: 1 = x, 2 = y, 0 = both
% 9/14/21

distance = 0;
for i = 1:(coord_idx2-coord_idx1)
    curr_i = coord_idx1 + i;
    prev_i = curr_i - 1;

    x_curr = mocap(curr_i,1);
    y_curr = mocap(curr_i,3);
    x_prev = mocap(prev_i,1);
    y_prev = mocap(prev_i,3);

    if axis == 1
        y_curr = 0;
        y_prev = 0;
    elseif axis == 2
        x_curr = 0;
        x_prev = 0;
    end

    new_distance = sqrt( (x_curr-x_prev)^2 + (y_curr-y_prev)^2 );
    if ~isnan(new_distance)
        distance = distance + new_distance;
    end
end

