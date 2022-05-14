function [estimateID, step_times] = recursive_stepID(curr_ID_labels, real_impact_times, step_times, step_o, step_x)
% recursively adds and takes out step times to match best step time outcome
% used by experiment4_footfalldetection.m
% 4/2/22

powerN = length(step_times); % unsemicolon to debug ################

allcombinations = ones(1, powerN); % no x's

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

for i = 1:allcombN
    curr_arr = allcombinations(i,:);
    one_idx = find(curr_arr == 1);
    two_idx = find(curr_arr == 2);
    score = 0;
    
    for j = 2:length(one_idx)
        curr_i = one_idx(j);
        past_i = one_idx(j-1);
        curr_diff = abs(step_times(curr_i) - step_times(past_i));
        diff_matrix(i,curr_i) = curr_diff;
        score = score + abs(curr_diff - step_o);
    end
    for k = 2:length(two_idx)
        curr_i = two_idx(k);
        past_i = two_idx(k-1);
        curr_diff = abs(step_times(curr_i) - step_times(past_i));
        diff_matrix(i,curr_i) = curr_diff;
        score = score + abs(curr_diff - step_x);
    end
    scores(i) = score;
end

[~,idx] = min(scores);
estimateID = allcombinations(idx,:);
% if first four impacts are all from the same person
while isempty(find(estimateID(1:4) == 1)) | isempty(find(estimateID(1:4) == 2))
    scores(idx) = max(scores); % ignore this combination
    [~,idx] = min(scores);
    estimateID = allcombinations(idx,:);
end

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
diff_matrix(idx,:) % uncomment to debug !!!!!!!!!!!!!!!!!
% diff_matrix(idx,est_o)
bigidxone = find(diff_matrix(idx,est_o) > step_o*1.45);
% diff_matrix(idx,est_x)
bigidxtwo = find(diff_matrix(idx,est_x) > step_x*1.45);
bigidx = sort([est_o(bigidxone), est_x(bigidxtwo)]);
if ~isempty(bigidx)
    [change_val,~] = max(bigidx);
    change_val % uncomment to debug!!!!!!!!
    % add element at bigidx-1 to duplicate it
    if step_times(bigidx(1)) == step_times(bigidx(1)-1) % if current one is overlap, add element at bigidx-2
        step_times = [step_times(1:bigidx(1)-2); step_times(bigidx(1)-2:end)];
    else
        step_times = [step_times(1:bigidx(1)-1); step_times(bigidx(1)-1:end)];
    end
    [estimateID, step_times] = recursive_stepID(curr_ID_labels, real_impact_times, step_times, step_o, step_x);
else
    % if all element in a column in diff_matrix is too small, delete, recursive
    % start from second element in both bc first element has diff = 0
    smallidxone = find(diff_matrix(idx,est_o(2:end)) < step_o/1.5);
    smallidxtwo = find(diff_matrix(idx,est_x(2:end)) < step_x/1.5);
    smallidx = sort([est_o(smallidxone+1), est_x(smallidxtwo+1)]); % weird that these are rows, bigidx is cols
    if ~isempty(smallidx)
        % delete element at smallidx
        step_times = [step_times(1:smallidx(1)-1); step_times(smallidx(1)+1:end)];
        [estimateID, step_times] = recursive_stepID(curr_ID_labels, real_impact_times, step_times, step_o, step_x);
    end
end


end

