function [est_impacts,iter] = overlapping_step_labeling_limp(est_impacts,detected_seg_starts,detected_seg_ends,step_o1,step_o2,step_x1,is_startBig,Fs_pcb)
% to find overlapping impacts and change labels
% 6/13/22 used in experiment4_processing3.m
% for limp 6/15/22

overlap_thresh = 0.36;
overlap_wrong = true;
o_idx = find(est_impacts(:,2) == 1);
x_idx = find(est_impacts(:,2) == 2);
iter = 0;
while overlap_wrong
    overlap_wrong = false;
    for i = 1:length(est_impacts(:,2))
        if ~isempty(find(round(detected_seg_starts.*Fs_pcb./10) == round(est_impacts(i,1)/10),1))
            % it's the start of the segment
            o_first_idx = find(est_impacts(o_idx,1) >= est_impacts(i,1),1);
            o_first = est_impacts(o_idx(o_first_idx),1);
            x_first_idx = find(est_impacts(x_idx,1) >= est_impacts(i,1),1);
            x_first = est_impacts(x_idx(x_first_idx),1);
            changed = false;
            if o_first < x_first
                % o happened first
                o_last_idx = find(est_impacts(o_idx,1) <= x_first & est_impacts(o_idx,1) >= o_first);
                % loop through all o's before x to see if overlapping x fit
                for j = 1:length(o_last_idx)
                    curr_idx = o_idx(o_last_idx(end-j + 1));
                    o_last = est_impacts(curr_idx,1);
                    if abs((x_first - o_last)/Fs_pcb - step_x1) < step_x1*overlap_thresh
                        % better to make overlapping
                        new_array = [est_impacts(curr_idx,1),2,est_impacts(curr_idx,3:end)];
                        est_impacts = [est_impacts(1:curr_idx,:);new_array;est_impacts(curr_idx+1:end,:)];
                        overlap_wrong = true;
                        changed = true;
                        o_idx = find(est_impacts(:,2) == 1);
                        x_idx = find(est_impacts(:,2) == 2);
                        iter = iter + 1;
                        break;
                    end
                end
            elseif o_first > x_first
                % x happened first
                x_last_idx = find(est_impacts(x_idx,1) <= o_first & est_impacts(x_idx,1) >= x_first);
                % loop through all o's before x to see if overlapping x fit
                for j = 1:length(x_last_idx)
                    curr_idx = x_idx(x_last_idx(end-j + 1));
                    x_last = est_impacts(curr_idx,1);
                    if is_startBig 
                        step_o = min(step_o1,step_o2);
                        if abs((o_first - x_last)/Fs_pcb - step_o) < step_o*overlap_thresh
                            % better to make overlapping
                            new_array = [est_impacts(curr_idx,1),1,est_impacts(curr_idx,3:end)];
                            est_impacts = [est_impacts(1:curr_idx,:);new_array;est_impacts(curr_idx+1:end,:)];
                            overlap_wrong = true;
                            changed = true;
                            o_idx = find(est_impacts(:,2) == 1);
                            x_idx = find(est_impacts(:,2) == 2);
                            iter = iter + 1;
                            is_startBig = false;
                            break;
                        end
                    else
                        step_o = min(step_o1,step_o2);
                        if abs((o_first - x_last)/Fs_pcb - step_o) < step_o*overlap_thresh
                            % better to make overlapping
                            new_array = [est_impacts(curr_idx,1),1,est_impacts(curr_idx,3:end)];
                            est_impacts = [est_impacts(1:curr_idx,:);new_array;est_impacts(curr_idx+1:end,:)];
                            overlap_wrong = true;
                            changed = true;
                            o_idx = find(est_impacts(:,2) == 1);
                            x_idx = find(est_impacts(:,2) == 2);
                            iter = iter + 1;
                            is_startBig = true;
                            break;
                        end
                    end
                end
            end
            if changed
                break;
            end
        end
        % check for ends of segments
        if ~isempty(find(round(detected_seg_ends.*Fs_pcb./10) == round(est_impacts(i,1)/10),1))
            % it's the end of the segment
            o_last_idx = find(est_impacts(o_idx,1) <= est_impacts(i,1));
            o_last = est_impacts(o_idx(o_last_idx(end)),1);
            x_last_idx = find(est_impacts(x_idx,1) <= est_impacts(i,1));
            x_last = est_impacts(x_idx(x_last_idx(end)),1);
            changed = false;
            if o_last > x_last
                % o happened last
                o_last_idx = find(est_impacts(o_idx,1) >= x_last & est_impacts(o_idx,1) <= o_last);
                % loop through all o's before x to see if overlapping x fit
                for j = 1:length(o_last_idx)
                    curr_idx = o_idx(o_last_idx(j));
                    curr_last = est_impacts(curr_idx,1);
                    if abs((curr_last - x_last)/Fs_pcb - step_x1) < step_x1*overlap_thresh
                        % better to make overlapping
                        new_array = [est_impacts(curr_idx,1),2,est_impacts(curr_idx,3:end)];
                        est_impacts = [est_impacts(1:curr_idx,:);new_array;est_impacts(curr_idx+1:end,:)];
                        overlap_wrong = true;
                        changed = true;
                        o_idx = find(est_impacts(:,2) == 1);
                        x_idx = find(est_impacts(:,2) == 2);
                        iter = iter + 1;
                        break;
                    end
                end
            elseif o_last < x_last
                % x happened last
                x_last_idx = find(est_impacts(x_idx,1) >= o_last & est_impacts(x_idx,1) <= x_last);
                % loop through all o's before x to see if overlapping x fit
                for j = 1:length(x_last_idx)
                    curr_idx = x_idx(x_last_idx(j));
                    curr_last = est_impacts(curr_idx,1);
                    if is_startBig
                        if mod(o_idx(o_last_idx(end))+1,2) == 0
                            step_o = max(step_o1,step_o2);
                        else
                            step_o = min(step_o1,step_o2);
                        end
                    else
                        if mod(o_idx(o_last_idx(end))+1,2) == 0
                            step_o = min(step_o1,step_o2);
                        else
                            step_o = max(step_o1,step_o2);
                        end
                    end
                    if abs((curr_last - o_last)/Fs_pcb - step_o) < step_o*overlap_thresh
                        % better to make overlapping
                        new_array = [est_impacts(curr_idx,1),1,est_impacts(curr_idx,3:end)];
                        est_impacts = [est_impacts(1:curr_idx,:);new_array;est_impacts(curr_idx+1:end,:)];
                        overlap_wrong = true;
                        changed = true;
                        o_idx = find(est_impacts(:,2) == 1);
                        x_idx = find(est_impacts(:,2) == 2);
                        iter = iter + 1;
                        break;
                    end
                end
            end
            if changed
                break;
            end
        end
    end
end

end
