function [estimateID_LR, step_times] = recursive_stepID_limp(curr_ID_labels,real_impact_times,step_times,step_o1,step_o2,step_x1,step_x2)
% recursively adds and takes out step times to match best step time outcome
% used by experiment4_footfalldetection.m &
% experiment4_allprocessingcompiled.m
% incorporates having two step times per person in case of limp
% 4/26/22

powerN = length(step_times); % unsemicolon to debug ################

allcombinations = ones(1, powerN); % no x's
label_allcombinations = ones(1,powerN); % same as allcombinations but left vs right on each person labeled
% 1 = x1, 2 = x2, 3 = o1, 4 = o2

for i = round(powerN/3):round(2*powerN/3) % know they at least have to be 33% of total footsteps
    curr_arr = ones(1,powerN);
    curr_arr(1:i) = 2;
    P = uniqueperms(curr_arr);
    allcombinations = [allcombinations; P];
end

% below uses the implemented dfs algorithm
% P1 = dfs_IDsequence(1,powerN,[1],step_times,step_times(1),step_times(2),1,0,zeros(1,powerN));
% P1(1,:) = []; % delete initialization
% allcombinations = [allcombinations; P1];
% P2 = dfs_IDsequence(1,powerN,[2],step_times,step_times(1),step_times(2),0,1,zeros(1,powerN));
% P2(1,:) = []; % delete initialization
% allcombinations = [allcombinations; P2];
% size(allcombinations)

allcombinations(1,:) = []; % delete first row from initialization
[allcombN,~] = size(allcombinations);
scores = zeros(allcombN,1);
diff_matrix = zeros(allcombN, powerN);

% accumulate scores for each possible combination
for i = 1:allcombN
    curr_arr = allcombinations(i,:);
    one_idx = find(curr_arr == 1);
    two_idx = find(curr_arr == 2);
    score = 0;
    
    % try starting with _o1 and _o2 to get 2 scores, choose better one
    score1 = 0;
    was_one = false;
    for j = 2:length(one_idx)
        curr_i = one_idx(j);
        past_i = one_idx(j-1);
        curr_diff = abs(step_times(curr_i) - step_times(past_i));
        diff_matrix(i,curr_i) = curr_diff;
        if was_one % used to alternate
            score1 = score1 + abs(curr_diff - step_o2);
            was_one = false;
        else
            score1 = score1 + abs(curr_diff - step_o1);
            was_one = true;
        end
    end
    
    score2 = 0;
    was_one = true; % this is the difference
    for j = 2:length(one_idx)
        curr_i = one_idx(j);
        curr_diff = diff_matrix(i,curr_i);
        if was_one % used to alternate
            score2 = score2 + abs(curr_diff - step_o2);
            was_one = false;
        else
            score2 = score2 + abs(curr_diff - step_o1);
            was_one = true;
        end
    end
    
    [minval, minidx_o] = min([score1, score2]);
    score = score + minval;
    
    % same thing for other person
    score1 = 0;
    was_one = false;
    for k = 2:length(two_idx)
        curr_i = two_idx(k);
        past_i = two_idx(k-1);
        curr_diff = abs(step_times(curr_i) - step_times(past_i));
        diff_matrix(i,curr_i) = curr_diff;
        if was_one % used to alternate
            score1 = score1 + abs(curr_diff - step_x2);
            was_one = false;
        else
            score1 = score1 + abs(curr_diff - step_x1);
            was_one = true;
        end
    end
    
    score2 = 0;
    was_one = true;
    for k = 2:length(two_idx)
        curr_i = two_idx(k);
        curr_diff = diff_matrix(i,curr_i);
        if was_one % used to alternate
            score2 = score2 + abs(curr_diff - step_x2);
            was_one = false;
        else
            score2 = score2 + abs(curr_diff - step_x1);
            was_one = true;
        end
    end
    
    [minval, minidx_x] = min([score1, score2]);
    score = score + minval;
    
    scores(i) = score;
    
    if minidx_o == 2 % start with o1
        one1_idx = one_idx(1:2:end);
        one2_idx = one_idx(2:2:end);
    else % start with o2
        one1_idx = one_idx(2:2:end);
        one2_idx = one_idx(1:2:end);
    end
    if minidx_x == 2 % start with x1
        two1_idx = two_idx(1:2:end);
        two2_idx = two_idx(2:2:end);
    else % start with x2
        two1_idx = two_idx(2:2:end);
        two2_idx = two_idx(1:2:end);
    end
    
    label_allcombinations(i,one1_idx) = 1;
    label_allcombinations(i,one2_idx) = 2;
    label_allcombinations(i,two1_idx) = 3;
    label_allcombinations(i,two2_idx) = 4;
end

[~,idx] = min(scores);
estimateID = allcombinations(idx,:);
% if first four impacts are all from the same person
while isempty(find(estimateID(1:4) == 1)) | isempty(find(estimateID(1:4) == 2))
    scores(idx) = max(scores); % ignore this combination
    [~,idx] = min(scores);
    estimateID = allcombinations(idx,:);
end

estimateID_LR = label_allcombinations(idx,:);

% display current state of step times
figure;
real_o = find(curr_ID_labels == 1);
real_x = find(curr_ID_labels == 2);
plot(real_impact_times(real_o),0,'bo')
hold on
plot(real_impact_times(real_x),0,'bx')
ylim([-1 1])

est_o = find(allcombinations(idx,:) == 1);
est_x = find(allcombinations(idx,:) == 2);

plot(step_times(est_o),0.5,'ro')
plot(step_times(est_x),0.5,'rx')
set(gcf,'Position',[100 100 500 200])

% perform adding impact first (more accurate?)
% if any step time is 1.5x larger than largest step time, missing an impact
diff_matrix(idx,:) % uncomment to debug !
% diff_matrix(idx,est_o)
est_o1 = find(label_allcombinations(idx,:) == 1);
est_o2 = find(label_allcombinations(idx,:) == 2);
est_x1 = find(label_allcombinations(idx,:) == 3);
est_x2 = find(label_allcombinations(idx,:) == 4);

bigidxone1 = find(diff_matrix(idx,est_o1) > step_o1*1.5);
bigidxone2 = find(diff_matrix(idx,est_o2) > step_o2*1.5);

% diff_matrix(idx,est_x)
bigidxtwo1 = find(diff_matrix(idx,est_x1) > step_x1*1.5);
bigidxtwo2 = find(diff_matrix(idx,est_x2) > step_x2*1.5);

bigidx = sort([est_o1(bigidxone1);est_o2(bigidxone2);est_x1(bigidxtwo1);est_x2(bigidxtwo2)]);
if ~isempty(bigidx)
    [bad_idx,~] = max(bigidx);
    bad_idx % uncomment to debug!!!!!!!!!!!!!!!!
    curr_idx = bigidx(1);
    
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
    
    
    
    step_times = [step_times(1:change_idx);step_times(change_idx:end)];
    [estimateID_LR, step_times] = recursive_stepID_limp(curr_ID_labels, real_impact_times, step_times, step_o1, step_o2, step_x1, step_x2);

else%     % below is old code:
%     [change_val,~] = max(bigidx);
%     change_val % uncomment to debug!!!!!!!!
%     % add element at bigidx-1 to duplicate it
%     if step_times(bigidx(1)) == step_times(bigidx(1)-1) % if current one is overlap, add element at bigidx-2
%         step_times = [step_times(1:bigidx(1)-2); step_times(bigidx(1)-2:end)];
%     else
%         step_times = [step_times(1:bigidx(1)-1); step_times(bigidx(1)-1:end)];
%     end
%     [estimateID_LR, step_times] = recursive_stepID_limp(curr_ID_labels,real_impact_times,step_times,step_o1,step_o2,step_x1,step_x2);

    % if all element in a column in diff_matrix is too small, delete, recursive
    % start from second element in both bc first element has diff = 0
    smallidxone1 = find(diff_matrix(idx,est_o1(2:end)) < step_o1/1.55);
    smallidxone2 = find(diff_matrix(idx,est_o2(2:end)) < step_o2/1.55);
    smallidxtwo1 = find(diff_matrix(idx,est_x1(2:end)) < step_x1/1.55);
    smallidxtwo2 = find(diff_matrix(idx,est_x2(2:end)) < step_x2/1.55);
    smallidx = sort([est_o1(smallidxone1+1),est_o2(smallidxone2+1),est_x1(smallidxtwo1+1),est_x2(smallidxtwo2+1)]); % weird that these are rows, bigidx is cols
    if ~isempty(smallidx)
        % delete element at smallidx
        step_times = [step_times(1:smallidx(1)-1); step_times(smallidx(1)+1:end)];
        [estimateID_LR, step_times] = recursive_stepID_limp(curr_ID_labels,real_impact_times,step_times,step_o1,step_o2,step_x1,step_x2);
    end
end


end

