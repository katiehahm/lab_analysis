function [diff_array, label_array, estimateID] = notrecursive_stepID_limp_whole(curr_ID_labels,real_impact_times,step_times,step_times_scores,step_o1,step_o2,step_x1,step_x2,plotbool)
% make sure step_times is 22 long max
% recursion happens in parent code using while loop
% adds and takes out step times to match best step time outcome
% used by experiment4_footfalldetection.m &
% experiment4_allprocessingcompiled.m
% incorporates having two step times per person in case of limp
% 5/19/22

powerN = length(step_times); % unsemicolon to debug ################

% allcombinations = ones(1, powerN); % no x's
% label_allcombinations = ones(1,powerN); % same as allcombinations but left vs right on each person labeled
% 1 = x1, 2 = x2, 3 = o1, 4 = o2

% calculate how big initialization should be
low_powerN = round(powerN/3); % know they at least have to be 33% of total footsteps
high_powerN = round(2*powerN/3);
comb_len = 0;
for i = low_powerN:high_powerN
    comb_len = comb_len + nchoosek(powerN,i);
end
allcombinations = ones(comb_len, powerN);
label_allcombinations = ones(comb_len, powerN);
% get all unique perms
assign_idx = 1;
for i = low_powerN:high_powerN 
    curr_arr = ones(1,powerN);
    curr_arr(1:i) = 2;
    P = uniqueperms(curr_arr);
    Prows = size(P,1);
    allcombinations(assign_idx:assign_idx+Prows-1,:) = P;
    assign_idx = assign_idx + Prows;
end

% below uses the implemented dfs algorithm
% P1 = dfs_IDsequence(1,powerN,[1],step_times,step_times(1),step_times(2),1,0,zeros(1,powerN));
% P1(1,:) = []; % delete initialization
% allcombinations = [allcombinations; P1];
% P2 = dfs_IDsequence(1,powerN,[2],step_times,step_times(1),step_times(2),0,1,zeros(1,powerN));
% P2(1,:) = []; % delete initialization
% allcombinations = [allcombinations; P2];
% size(allcombinations)

% allcombinations(1,:) = []; % delete first row from initialization
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
while ~any(estimateID(1:4) == 1 | estimateID(1:4) == 2)
% while isempty(find(estimateID(1:4) == 1)) || isempty(find(estimateID(1:4) == 2))
    scores(idx) = max(scores); % ignore this combination
    [~,idx] = min(scores);
    estimateID = allcombinations(idx,:);
end

diff_array = diff_matrix(idx,:);
label_array = label_allcombinations(idx,:);

if plotbool
    % display current state of step times
    figure;
    % real_o = find(curr_ID_labels == 1);
    % real_x = find(curr_ID_labels == 2);
    plot(real_impact_times(curr_ID_labels == 1),0,'bo')
    hold on
    plot(real_impact_times(curr_ID_labels == 2),0,'bx')
    ylim([-1 1])

    % est_o = find(allcombinations(idx,:) == 1);
    % est_x = find(allcombinations(idx,:) == 2);

    plot(step_times(allcombinations(idx,:) == 1),0.5,'ro')
    plot(step_times(allcombinations(idx,:) == 2),0.5,'rx')
    set(gcf,'Position',[100 100 500 200])
end

end

