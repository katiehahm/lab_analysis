function [est_impacts,fixed] = big_step_labeling_x(est_impacts,detected_seg_starts,avg_step, avg_other_step,Fs_pcb)
% this is to clean up code in experiment4_processing3.m
% to add or replace person labeling if there's a big step
% 6/10/22

o_idx = find(est_impacts(:,2) == 1);
x_idx = find(est_impacts(:,2) == 2);
fixed = false;
for i = 2:length(x_idx)
    curr_time = est_impacts(x_idx(i),1)/Fs_pcb;
    prev_time = est_impacts(x_idx(i-1),1)/Fs_pcb;
    if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_time,1) % in same segment
        if curr_time - prev_time > 1.75*avg_step
            fixed = true;
            times_inbetween = est_impacts(x_idx(i-1):x_idx(i),1)./Fs_pcb - prev_time;
            [~,minidx] = min(abs(times_inbetween - avg_step));
            % check if this should be replaced with x or made overlapping
            curr_idx = x_idx(i-1) + minidx - 1;
            curr_o_idx = find(o_idx == curr_idx);
            next_o_time = est_impacts(o_idx(curr_o_idx)+1,1)./Fs_pcb;
            prev_o_time = est_impacts(o_idx(curr_o_idx)-1,1)./Fs_pcb;
            curr_o_time = est_impacts(o_idx(curr_o_idx),1)./Fs_pcb;
            if abs((next_o_time-prev_o_time) - avg_other_step) < abs((next_o_time-curr_o_time) - avg_other_step)
                % replace
                est_impacts(curr_idx,2) = 2;
            else
                % overlapping
                new_array = [est_impacts(curr_idx,1),2,est_impacts(curr_idx,3:end)];
                est_impacts = [est_impacts(1:curr_idx,:);new_array;est_impacts(curr_idx+1:end,:)];
            end
            break;
        end
    end
end
end

