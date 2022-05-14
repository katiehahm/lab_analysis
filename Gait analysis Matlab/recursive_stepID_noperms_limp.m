function [estimateID_LR, step_times] = recursive_stepID_noperms_limp(curr_ID_labels, real_impact_times, estimateID_LR, step_times, step_o1, step_o2, step_x1, step_x2)
% recursively adds and takes out step times to match best step time outcome
% same as recursive_stepID_noperms, but includes limps

% perform adding impact first (more accurate?)
% if any step time is 1.5x larger than largest step time, missing an impact
est_o = find(estimateID_LR == 1 | estimateID_LR == 2);
est_x = find(estimateID_LR == 3 | estimateID_LR == 4);
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

diff_matrix % uncomment to debug!!!!!!!!!!!!!!!!!!!!!!

% redo L/R labeling bc adding/removing element changes it up
score = zeros(1,4);
o_odd_idx = est_o(1:2:end);
o_eve_idx = est_o(2:2:end);
x_odd_idx = est_x(1:2:end);
x_eve_idx = est_x(2:2:end);

o_odd_start = sum(abs(diff_matrix(o_odd_idx)-step_o1)) + sum(abs(diff_matrix(o_eve_idx)-step_o2));
o_eve_start = sum(abs(diff_matrix(o_odd_idx)-step_o2)) + sum(abs(diff_matrix(o_eve_idx)-step_o1));
x_odd_start = sum(abs(diff_matrix(x_odd_idx)-step_x1)) + sum(abs(diff_matrix(x_eve_idx)-step_x2));
x_eve_start = sum(abs(diff_matrix(x_odd_idx)-step_x2)) + sum(abs(diff_matrix(x_eve_idx)-step_x1));

if o_odd_start < o_eve_start
    estimateID_LR(o_odd_idx) = 1;
    estimateID_LR(o_eve_idx) = 2;
else
    estimateID_LR(o_odd_idx) = 2;
    estimateID_LR(o_eve_idx) = 1;
end
if x_odd_start < x_eve_start
    estimateID_LR(x_odd_idx) = 3;
    estimateID_LR(x_eve_idx) = 4;
else
    estimateID_LR(x_odd_idx) = 4;
    estimateID_LR(x_eve_idx) = 3;
end


% display current state of step times
figure;
plot(real_impact_times(curr_ID_labels == 1),0,'bo')
hold on
plot(real_impact_times(curr_ID_labels == 2),0,'bx')
ylim([-1 2])
hold on
est_o1 = find(estimateID_LR == 1);
est_o2 = find(estimateID_LR == 2);
est_x1 = find(estimateID_LR == 3);
est_x2 = find(estimateID_LR == 4);
plot(step_times(est_o1),1,'ro')
plot(step_times(est_o2),1,'ro')
plot(step_times(est_x1),1,'rx')
plot(step_times(est_x2),1,'rx')
set(gcf,'Position',[100 100 500 200])



bigidxone1 = find(diff_matrix(est_o1) > step_o1*1.55);
bigidxone2 = find(diff_matrix(est_o2) > step_o2*1.55);
bigidxtwo1 = find(diff_matrix(est_x1) > step_x1*1.55);
bigidxtwo2 = find(diff_matrix(est_x2) > step_x2*1.55);
bigidx = sort([est_o1(bigidxone1);est_o2(bigidxone2);est_x1(bigidxtwo1);est_x2(bigidxtwo2)]);
if ~isempty(bigidx)
    [bad_idx,~] = max(bigidx);
    bad_idx % uncomment to debug!!!!!!!!!!!!!!!!
%     curr_idx = bigidx(1);
    curr_idx = bad_idx;
    
    if estimateID_LR(curr_idx) == 1
        [~,change_idx] = min(abs(step_times(curr_idx)-step_times(1:curr_idx-1) - step_o1));
    elseif estimateID_LR(curr_idx) == 2
        [~,change_idx] = min(abs(step_times(curr_idx)-step_times(1:curr_idx-1) - step_o2));
    elseif estimateID_LR(curr_idx) == 3
        [~,change_idx] = min(abs(step_times(curr_idx)-step_times(1:curr_idx-1) - step_x1));
    else
        [~,change_idx] = min(abs(step_times(curr_idx)-step_times(1:curr_idx-1) - step_x2));
    end
    change_idx
    
%     if estimateID_LR(curr_idx) == 1 | estimateID_LR(curr_idx) == 2
%         past_idx = curr_idx -1;
%         while(estimateID_LR(past_idx) == 1 | estimateID_LR(past_idx) == 2)
%             past_idx = past_idx - 1;
%         end
%         if estimateID_LR(past_idx) == 3
%             estimateID = [estimateID_LR(1:change_idx),4,estimateID_LR(change_idx+1:end)];
%         else
%             estimateID = [estimateID_LR(1:change_idx),3,estimateID_LR(change_idx+1:end)];
%         end
%     else
%         past_idx = change_idx -1;
%         while(estimateID_LR(past_idx) == 3 | estimateID_LR(past_idx) == 4)
%             past_idx = past_idx - 1;
%         end
%         if estimateID_LR(past_idx) == 1
%             estimateID = [estimateID_LR(1:change_idx),2,estimateID_LR(change_idx+1:end)];
%         else
%             estimateID = [estimateID_LR(1:change_idx),1,estimateID_LR(change_idx+1:end)];
%         end
%     end
    if estimateID_LR(curr_idx) == 1
        new_ID = 2;
    elseif estimateID_LR(curr_idx) == 2
        new_ID = 1;
    elseif estimateID_LR(curr_idx) == 3
        new_ID = 4;
    else
        new_ID = 3;
    end
    
    estimateID_LR = [estimateID_LR(1:change_idx),new_ID,estimateID_LR(change_idx+1:end)];
    step_times = [step_times(1:change_idx);step_times(change_idx:end)];
    [estimateID_LR, step_times] = recursive_stepID_noperms_limp(curr_ID_labels, real_impact_times, estimateID_LR, step_times, step_o1, step_o2, step_x1, step_x2);

else
    % if all element in a column in diff_matrix is too small, delete, recursive
    % start from second element in both bc first element has diff = 0
    smallidxone1 = find(diff_matrix(est_o1(2:end)) < step_o1/1.55); % THIS IS NOT 2
    smallidxone2 = find(diff_matrix(est_o2(2:end)) < step_o2/1.55);
    smallidxtwo1 = find(diff_matrix(est_x1(2:end)) < step_x1/1.55);
    smallidxtwo2 = find(diff_matrix(est_x2(2:end)) < step_x2/1.55);
    smallidx = sort([est_o1(smallidxone1+1),est_o2(smallidxone2+1),est_x1(smallidxtwo1+1),est_x2(smallidxtwo2+1)]); % weird that these are rows, bigidx is cols
    if ~isempty(smallidx)
        % delete element at smallidx
        step_times = [step_times(1:smallidx(1)-1); step_times(smallidx(1)+1:end)];
        estimateID_LR = [estimateID_LR(1:smallidx(1)-1), estimateID_LR(smallidx(1)+1:end)];
        [estimateID_LR, step_times] = recursive_stepID_noperms_limp(curr_ID_labels, real_impact_times, estimateID_LR, step_times, step_o1, step_o2, step_x1, step_x2);
    end
end


end

