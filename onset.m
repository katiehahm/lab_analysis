function [disturb, disturb_val] = onset(p_idx, filt)
% Given the peaks of each hit, finds the onset of that disturbance
% 03/20/20
npeaks = length(p_idx);
disturb = zeros(npeaks, 4);
disturb_val = zeros(npeaks, 4);
tolerance = 5000; % max # samples between onset and peak
tol_val = 0.001; % min value to surpass to qualify as signal
for c = 1:4
    for i = 1:npeaks
        for j = 1:tolerance
            curr_i = p_idx(i, c) - tolerance + j;
            if curr_i < 1
                curr_i = 1;
            end
            val = filt(curr_i, c);
            if abs(val) > tol_val
                disturb(i,c) = curr_i;
                break;
            end
        end
    end
end

for c = 1:4
    for i = 1:npeaks
        disturb_val(i,c) = filt(disturb(i,c),c);
    end
end

end

