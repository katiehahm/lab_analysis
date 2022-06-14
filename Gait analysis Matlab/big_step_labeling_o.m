function [est_impacts,fixed] = big_step_labeling_o(est_impacts,detected_seg_starts,avg_step,curr_person_label,Fs_pcb)
% this is to clean up code in experiment4_processing3.m
% to add or replace person labeling if there's a big step person o
% 6/10/22

% if curr_person_label = 1, then it's person o
% if " = 2, then it's person x

fixed = false;
idx = find(est_impacts(:,2) == curr_person_label);
for i = 2:length(idx)
    curr_time = est_impacts(idx(i),1)/Fs_pcb;
    prev_time = est_impacts(idx(i-1),1)/Fs_pcb;
    if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_time,1) % in same segment
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