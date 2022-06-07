function [diff_array, label_array, estimateID] = optimized_stepID(curr_ID_labels,real_impact_times,step_times,step_o1,step_o2,step_x1,step_x2)
% DO NOT use. This actually takes 2x long to run compared to notrecursive_stepID_limp_whole
% this is variation of notrecursive_stepID_limp_whole.m, but
% trying to use matrix alg to shorten runtime
% make sure step_times is 22 long max
% recursion happens in parent code using while loop
% adds and takes out step times to match best step time outcome
% used by experiment4_footfalldetection.m &
% experiment4_allprocessingcompiled.m
% incorporates having two step times per person in case of limp
% 5/20/22

powerN = length(step_times); % unsemicolon to debug ################

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

% allcombinations(1,:) = []; % delete first row from initialization
% [allcombN,~] = size(allcombinations);
scores = zeros(comb_len,1);
person1_idx = allcombinations == 1;
person2_idx = allcombinations == 2;
diff_matrix = zeros(comb_len, powerN);

% accumulate scores for each row of possible combination
for i = 1:comb_len
    one_arr = step_times(person1_idx(i,:)).';
    one_diff = [0,diff(one_arr)]; % insert 0 to keep array lengths same
    diff_matrix(i,person1_idx(i,:)) = one_diff;
    o1_odd = sum(abs(one_diff(1:2:end) - step_o1)) + sum(abs(one_diff(2:2:end) - step_o2));
    o1_even = sum(abs(one_diff(1:2:end) - step_o2)) + sum(abs(one_diff(2:2:end) - step_o1));
    [one_minval, one_minidx] = min([o1_odd,o1_even]);
    
    two_arr = step_times(person2_idx(i,:)).';
    two_diff = [0,diff(two_arr)];
    diff_matrix(i,person2_idx(i,:)) = two_diff;
    x1_odd = sum(abs(two_diff(1:2:end) - step_x1)) + sum(abs(two_diff(2:2:end) - step_x2));
    x1_even = sum(abs(two_diff(1:2:end) - step_x2)) + sum(abs(two_diff(2:2:end) - step_x1));
    [two_minval, two_minidx] = min([x1_odd,x1_even]);
    
    scores(i) = one_minval + two_minval;
    
    one_idx = find(person1_idx(i,:) == 1);
    one_odd_idx = one_idx(1:2:end);
    one_even_idx = one_idx(2:2:end);
    two_idx = find(person2_idx(i,:) == 1);
    two_odd_idx = two_idx(1:2:end);
    two_even_idx = two_idx(2:2:end);
    
    label_allcombinations(i,one_odd_idx) = one_minidx;
    label_allcombinations(i,one_even_idx) = floor((one_minidx+1)/3) + mod(one_minidx+1,3); % = 2 if one_minidx = 1, =1 if " = 2
    label_allcombinations(i,two_odd_idx) = 2 + two_minidx;
    label_allcombinations(i,two_even_idx) = 2 + floor((two_minidx+1)/3) + mod(two_minidx+1,3);
end

[~,idx] = min(scores);
diff_array = diff_matrix(idx,:);
label_array = label_allcombinations(idx,:);


estimateID = allcombinations(idx,:);
% if first four impacts are all from the same person
while ~any(estimateID(1:4) == 1 | estimateID(1:4) == 2)
% while isempty(find(estimateID(1:4) == 1)) || isempty(find(estimateID(1:4) == 2))
    scores(idx) = max(scores); % ignore this combination
    [~,idx] = min(scores);
    estimateID = allcombinations(idx,:);
end

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