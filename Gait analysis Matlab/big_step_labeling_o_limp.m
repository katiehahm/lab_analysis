function [est_impacts,fixed] = big_step_labeling_o_limp(est_impacts,detected_seg_starts,step_o1,step_o2,Fs_pcb,startBig)
% this is to clean up code in experiment4_processing3.m
% to add or replace person labeling if there's a big step person o
% incorporates the possible limp
% startBig = true if first step was big step, false is first step was small
% 6/15/22


fixed = false;
idx = find(est_impacts(:,2) == 1);
for i = 2:length(idx)
    curr_time = est_impacts(idx(i),1)/Fs_pcb;
    prev_time = est_impacts(idx(i-1),1)/Fs_pcb;
    if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_time,1) % in same segment
        if mod(i,2) == 0 % is even
            if startBig
                avg_step = max(step_o1,step_o2);
            else
                avg_step = min(step_o1,step_o2);
            end
        else
            if startBig
                avg_step = min(step_o1,step_o2);
            else
                avg_step = max(step_o1,step_o2);
            end
        end
        if curr_time - prev_time > 1.75*avg_step
            fixed = true;
            times_inbetween = est_impacts(idx(i-1):idx(i),1)./Fs_pcb - prev_time;
            [~,minidx] = min(abs(times_inbetween - avg_step));
            fix_idx = idx(i-1) + minidx - 1;
            est_impacts(fix_idx,2) = 1;
            break;
        end
    end
end

end