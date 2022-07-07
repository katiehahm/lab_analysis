function [est_impacts,fixed] = small_step_labeling(est_impacts,detected_seg_starts,step_curr,step_other,curr_person,other_person,Fs_pcb)
% 6/13/22
% delete any small steps 

small_scaler = 0.4;

curr_idx = find(est_impacts(:,2) == curr_person);
other_idx = find(est_impacts(:,2) == other_person);
fixed = false;
for i = 2:length(curr_idx)
    curr_time = est_impacts(curr_idx(i),1)/Fs_pcb;
    prev_time = est_impacts(curr_idx(i-1),1)/Fs_pcb;
    if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_time,1) % in same segment
        if curr_time - prev_time < small_scaler*step_curr % found small step
            fixed = true;
            % decide if remove or relabel for both impacts
            x_prev = find(est_impacts(other_idx,1) <= curr_time*Fs_pcb);
            if ~isempty(x_prev) 
                x_prev = x_prev(end);
                x_prev_time = est_impacts(other_idx(x_prev),1)/Fs_pcb;
                if find(detected_seg_starts > x_prev_time,1) ~= find(detected_seg_starts > curr_time,1) % not in same seg
                    x_prev_time = 0;
                end
            else % prev does not exist
                x_prev_time = 0;
            end
            x_next = find(est_impacts(other_idx,1) >= curr_time*Fs_pcb);
            if ~isempty(x_next)
                x_next = x_next(1);
                x_next_time = est_impacts(other_idx(x_next),1)/Fs_pcb;
                if find(detected_seg_starts > x_next_time,1) ~= find(detected_seg_starts > curr_time,1) % not in same seg
                    x_next_time = 0;
                end
            else % next does not exist
                x_next_time = 0;
            end
            
            if x_prev_time == 0 % prev time doesn't exist
                curr_margin = abs(x_next_time - curr_time - step_other);
                prev_margin = abs(x_next_time - prev_time - step_other);
                if curr_margin < prev_margin
                    if curr_margin < step_other*0.3 
                        % better to replace current impact
                        curr_or_prev = 1;
                    else
                        % better to remove impact
                        curr_or_prev = 0;
                    end
                else
                    if prev_margin < step_other*0.3
                        % better to replace prev impact
                        curr_or_prev = 2;
                    else
                        % better to remove impact
                        curr_or_prev = 0;
                    end
                end
            elseif x_next_time == 0 % next time doesn't exist
                curr_margin = abs(curr_time - x_prev_time - step_other);
                prev_margin = abs(prev_time - x_prev_time - step_other);
                if curr_margin < prev_margin
                    if curr_margin < step_other*0.3 
                        % better to replace current impact
                        curr_or_prev = 1;
                    else
                        % better to remove impact
                        curr_or_prev = 0;
                    end
                else
                    if prev_margin < step_other*0.3
                        % better to replace prev impact
                        curr_or_prev = 2;
                    else
                        % better to remove impact
                        curr_or_prev = 0;
                    end
                end
            else
                % both prev and next time exist
                curr_margin = abs(x_next_time - curr_time - step_other) + abs(curr_time - x_prev_time - step_other);
                prev_margin = abs(x_next_time - prev_time - step_other) + abs(prev_time - x_prev_time - step_other);
                keep_margin = abs(x_next_time - x_prev_time - step_other);
                if curr_margin < keep_margin
                    % better to replace current impact
                    curr_or_prev = 1;
                elseif prev_margin < keep_margin
                    % better to replace prev impact
                    curr_or_prev = 2;
                else
                    % better to remove impact
                    curr_or_prev = 0;
                end
            end
                
                
            % check that prev and next x are in same seg
            
%             if find(detected_seg_starts > x_next_time,1) == find(detected_seg_starts > curr_time,1) % if in same segment
%                 if find(detected_seg_starts > x_prev_time,1) == find(detected_seg_starts > curr_time,1) % if in same segment
%                     curr_or_prev = 0;
%                     curr_margin = abs(x_next_time - curr_time - step_other) + abs(curr_time - x_prev_time - step_other);
%                     prev_margin = abs(x_next_time - prev_time - step_other) + abs(prev_time - x_prev_time - step_other);
%                     % check if better to replace
%                     if abs(x_next_time - x_prev_time - step_other) > abs(x_next_time - curr_time - step_other)
%                         if abs(x_next_time - x_prev_time - step_other) > abs(curr_time - x_prev_time - step_other)
%                             curr_or_prev = 1; % better to replace current impact
%                         end
%                     end
%                     if abs(x_next_time - x_prev_time - step_other) > abs(x_next_time - prev_time - step_other)
%                         if abs(x_next_time - x_prev_time - step_other) > abs(prev_time - x_prev_time - step_other)
%                             if curr_or_prev == 1 & prev_margin < curr_margin
%                                 curr_or_prev = 2; % better to replace previous impact
%                             elseif curr_or_prev == 0
%                                 curr_or_prev = 2;
%                             end
%                         end
%                     end
%                 else % only next x is in same seg
%                     curr_or_prev = 0;
%                     replaced_curr_time = abs(x_next_time - curr_time - step_other);
%                     replaced_prev_time = abs(x_next_time - prev_time - step_other);
%                     if replaced_prev_time < replaced_curr_time
%                         if replaced_prev_time < step_other*0.3 
%                             curr_or_prev = 2;
%                         end
%                     else
%                         if replaced_curr_time < step_other*0.3 
%                             curr_or_prev = 1;
%                         end
%                     end
%                 end
%             elseif find(detected_seg_starts > x_prev_time,1) == find(detected_seg_starts > curr_time,1)
%                 % only prev x is in same seg
%                 curr_or_prev = 0;
%                 replaced_curr_time = abs(curr_time - x_prev_time - step_other);
%                 replaced_prev_time = abs(prev_time - x_prev_time - step_other);
%                 if replaced_prev_time < replaced_curr_time
%                     if replaced_prev_time < step_other*0.3 
%                         curr_or_prev = 2;
%                     end
%                 else
%                     if replaced_curr_time < step_other*0.3 
%                         curr_or_prev = 1;
%                     end
%                 end
%             end
            % now edit est_impacts
            if curr_or_prev == 0 % remove impact
                if i < length(curr_idx)
                    next_o_time = est_impacts(curr_idx(i+1),1)/Fs_pcb;
                else
                    next_o_time = 0;
                end
                if i > 3
                    prev_o_time = est_impacts(curr_idx(i-2),1)/Fs_pcb;
                else
                    prev_o_time = 0;
                end
                if next_o_time == 0
                    if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_o_time,1)
                        curr_o_margin = abs(curr_time - prev_o_time - step_curr);
                        prev_o_margin = abs(prev_time - prev_o_time - step_curr);
                    end
                elseif prev_o_time == 0
                    if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > next_o_time,1) % in same segment
                        curr_o_margin = abs(next_o_time - curr_time - step_curr);
                        prev_o_margin = abs(next_o_time - prev_time - step_curr);
                    end
                else
                    if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > next_o_time,1) % in same segment
                        if find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_o_time,1) % in same segment
                            % both prev and next current person impact exists
                            curr_o_margin = abs(next_o_time - curr_time - step_curr) + abs(curr_time - prev_o_time - step_curr);
                            prev_o_margin = abs(next_o_time - prev_time - step_curr) + abs(prev_time - prev_o_time - step_curr);

                        else
                            % only next impact current person exists
                            curr_o_margin = abs(next_o_time - curr_time - step_curr);
                            prev_o_margin = abs(next_o_time - prev_time - step_curr);
                        end
                    elseif find(detected_seg_starts > curr_time,1) == find(detected_seg_starts > prev_o_time,1) % in same segment
                        % only prev impact current person exists
                        curr_o_margin = abs(curr_time - prev_o_time - step_curr);
                        prev_o_margin = abs(prev_time - prev_o_time - step_curr);
                    end
                end
                if curr_o_margin > prev_o_margin % delete current impact
                    est_impacts(curr_idx(i),:) = [];
                else % delete prev impact
                    est_impacts(curr_idx(i-1),:) = [];
                end
            elseif curr_or_prev == 1 % replace current impact
                est_impacts(curr_idx(i),2) = other_person;
            else % replace previous impact
                est_impacts(curr_idx(i-1),2) = other_person;
            end
            break;
        end
    end
end



end
