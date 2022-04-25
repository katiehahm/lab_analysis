function [estimateID, step_times] = recursive_stepID_noperms(curr_ID_labels, real_impact_times, estimateID, step_times, step_o, step_x)
% recursively adds and takes out step times to match best step time outcome
% to use after using recursive_stepID to find the best step times for two
% halves of segment. This is to make sure the middle segments are not empty
% because the ends of each segment are missing impacts
% used by experiment4_footfalldetection.m
% 4/23/22

% display current state of step times
figure;
real_o = find(curr_ID_labels == 1);
real_x = find(curr_ID_labels == 2);
plot(real_impact_times(real_o),0,'bo')
hold on
plot(real_impact_times(real_x),0,'bx')
ylim([-1 2])
hold on
est_o = find(estimateID == 1);
est_x = find(estimateID == 2);
plot(step_times(est_o),1,'ro')
plot(step_times(est_x),1,'rx')
set(gcf,'Position',[100 100 500 200])

% perform adding impact first (more accurate?)
% if any step time is 1.5x larger than largest step time, missing an impact

diff_matrix = zeros(1,length(step_times));
for i = 2:length(est_o)
    curr_i = est_o(i);
    past_i = est_o(i-1);
    curr_diff = abs(step_times(curr_i) - step_times(past_i));
    diff_matrix(curr_i) = curr_diff;
end
for i = 2:length(est_x)
    curr_i = est_x(i);
    past_i = est_x(i-1);
    curr_diff = abs(step_times(curr_i) - step_times(past_i));
    diff_matrix(curr_i) = curr_diff;
end

bigidxone = find(diff_matrix(est_o) > step_o*1.45);
bigidxtwo = find(diff_matrix(est_x) > step_x*1.45);
bigidx = sort([est_o(bigidxone); est_x(bigidxtwo)]);
if ~isempty(bigidx)
    [bad_idx,~] = max(bigidx);
%     bad_idx % uncomment to debug!!!!!!!!!!!!!!!!
    curr_idx = bigidx(1);
    
    if estimateID(curr_idx) == 1
        [~,change_idx] = min(abs(step_times(curr_idx)-step_times(1:curr_idx-1) - step_o));
    else
        [~,change_idx] = min(abs(step_times(curr_idx)-step_times(1:curr_idx-1) - step_x));
    end
    if estimateID(change_idx) == 1
        estimateID = [estimateID(1:change_idx),2,estimateID(change_idx+1:end)];
    else
        estimateID = [estimateID(1:change_idx),1,estimateID(change_idx+1:end)];
    end
    step_times = [step_times(1:change_idx);step_times(change_idx:end)];
    [estimateID, step_times] = recursive_stepID_noperms(curr_ID_labels, real_impact_times, estimateID, step_times, step_o, step_x);

else
    % if all element in a column in diff_matrix is too small, delete, recursive
    % start from second element in both bc first element has diff = 0
    smallidxone = find(diff_matrix(est_o(2:end)) < step_o/1.5);
    smallidxtwo = find(diff_matrix(est_x(2:end)) < step_x/1.5);
    smallidx = sort([est_o(smallidxone+1), est_x(smallidxtwo+1)]); % weird that these are rows, bigidx is cols
    if ~isempty(smallidx)
        % delete element at smallidx
        step_times = [step_times(1:smallidx(1)-1); step_times(smallidx(1)+1:end)];
        estimateID = [estimateID(1:smallidx(1)-1), estimateID(smallidx(1)+1:end)];
        [estimateID, step_times] = recursive_stepID_noperms(curr_ID_labels, real_impact_times, estimateID, step_times, step_o, step_x);
    end
end


end

