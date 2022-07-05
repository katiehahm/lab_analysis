% makes overall csv for alltakes localization, TA estimation, kmeans
% 6/10/22
%% getting rmse on GMM results

takes = {'regular1', 'limp1', 'limp2', 'weight1', 'weight2', 'regular2'};
exp_subject = 'Praneeth 5';
gmm_rmse = 0;
for t = 1:length(takes)
    filename = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\',exp_subject,'\ProcessedData\both_', char(takes(t)),'.mat'];
    load(filename)
    gmm_rmse = gmm_rmse + (estimated_scaled_means(1,1)-real_means(1,1))^2;
    gmm_rmse = gmm_rmse + (estimated_scaled_means(1,2)-real_means(1,2))^2;
    gmm_rmse = gmm_rmse + (estimated_scaled_means(2,1)-real_means(2,1))^2;
    gmm_rmse = gmm_rmse + (estimated_scaled_means(2,2)-real_means(2,2))^2;
    
end
gmm_rmse = sqrt(gmm_rmse/(4*length(takes)))

%% making a localization csv for all takes

takes = {'regular1', 'limp1', 'limp2', 'weight1', 'weight2', 'regular2'};
newA = zeros(1,32);
person = '2';
for t = 1:length(takes)
    filename = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\',exp_subject,'\ProcessedData\ExcelData\both_', char(takes(t)),'_localization_p',person','_withta.csv'];
    T = readtable(filename);
    A = table2array(T);
    newA = [newA;A];
end
newA(1,:) = []; % initialization

filename = ['C:/Users/Katie/Dropbox (MIT)/Lab/Analysis/Experiment4/',exp_subject,'/ProcessedData/ExcelData/alltakes_localization_p',person','_withta.csv'];
writematrix(newA,filename)

% then run exp4_notrecursive_localization.py for all takes

%% make TA csv with overlapping impacts 5/26/22 using alltakes localization results

filename1 = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\',exp_subject,'\ProcessedData\ExcelData\alltakes_localization_p1_results.csv'];
T = readtable(filename1);
loc_results1 = table2array(T);
filename2 = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\',exp_subject,'\ProcessedData\ExcelData\alltakes_localization_p2_results.csv'];
T = readtable(filename2);
loc_results2 = table2array(T);

takes = {'regular1', 'limp1', 'limp2', 'weight1', 'weight2', 'regular2'};
s1 = [-3.590,-3.343];
s2 = [-3.580,2.61];
s3 = [3.639,2.11];
s4 = [3.650,-3.412];
% these are values from interpolating s1-4
s5 = [3.61,2.36];
s6 = [3.62,-3.3775];

for p = 1:2
    person = int2str(p);
    rowcount = 1;
    for t = 1:length(takes)
        filename = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\',exp_subject,'\ProcessedData\ExcelData\both_', char(takes(t)),'_localization_p',person','_withta.csv'];
        T = readtable(filename);
        A = table2array(T);
        numrows = length(A(:,1));
        if p == 1
            pred_locs = loc_results1(:,2); % second col has est values
        else
            pred_locs = loc_results2(:,2);
        end
        locs = pred_locs(rowcount:rowcount + numrows - 1);
        featureV = zeros(1,37);
        for i = 1:numrows
            xcoord = locs(i);
            dist1 = sqrt( (xcoord-s1(1)).^2 );
            dist2 = sqrt( (xcoord-s2(1)).^2 );
            dist3 = sqrt( (xcoord-s3(1)).^2 );
            dist4 = sqrt( (xcoord-s4(1)).^2 );
            dist5 = sqrt( (xcoord-s5(1)).^2 );
            dist6 = sqrt( (xcoord-s6(1)).^2 );
            feature = [A(i,2),dist1,dist2,dist3,dist4,dist5,dist6,A(i,3:end)];
            featureV(end+1,:) = feature;
        end
        featureV(1,:) = []; % initialization
        newfilename = ['C:/Users/Katie/Dropbox (MIT)/Lab/Analysis/Experiment4/',exp_subject,'/ProcessedData/ExcelData/both_', char(takes(t)),'_ta_p',person,'.csv'];
        writematrix(featureV,newfilename)
        rowcount = rowcount + numrows;
    end
end

% then run exp4_TAestimation.py on each take

%% doing kmeans on TA estimation results all takes 5/24/22

takes = {'regular1', 'limp1', 'limp2', 'weight1', 'weight2', 'regular2'};
TA_rmse = 0;
figure;
hold on
for t = 1:length(takes)
    filename1 = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\',exp_subject,'\ProcessedData\ExcelData\both_', char(takes(t)),'_ta_p1_results.csv'];
    T1 = readtable(filename1);
    A1 = table2array(T1);
    est_TA1 = A1(2:end,3);
    [~,cent1] = kmeans(est_TA1,2);
    
    filename2 = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\',exp_subject,'\ProcessedData\ExcelData\both_', char(takes(t)),'_ta_p2_results.csv'];
    T2 = readtable(filename2);
    A2 = table2array(T2);
    est_TA2 = A2(2:end,3);
    [~,cent2] = kmeans(est_TA2,2);
    
    load(['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\',exp_subject,'\ProcessedData\both_',char(takes(t))])
    p1f1_idx = find(clean_impacts(:,2) == 11);
    p1f2_idx = find(clean_impacts(:,2) == 12);
    p2f1_idx = find(clean_impacts(:,2) == 21);
    p2f2_idx = find(clean_impacts(:,2) == 22);
    real_p1f1 = mean(clean_impacts(p1f1_idx,5));
    real_p1f2 = mean(clean_impacts(p1f2_idx,5));
    real_p2f1 = mean(clean_impacts(p2f1_idx,5));
    real_p2f2 = mean(clean_impacts(p2f2_idx,5));
    
    plot(min(cent1),min([real_p1f1,real_p1f2]),'ro')
    plot(max(cent1),max([real_p1f1,real_p1f2]),'bo')
    plot(min(cent2),min([real_p2f1,real_p2f2]),'rx')
    plot(max(cent2),max([real_p2f1,real_p2f2]),'bx')
    legend('Leg 1 person 1','Leg 2 person 1','Leg 1 person 2','Leg 2 person 2')
    
    % calculate rmse
    TA_rmse = TA_rmse + abs(min(cent1)-min([real_p1f1,real_p1f2]));
    TA_rmse = TA_rmse + abs(max(cent1)-max([real_p1f1,real_p1f2]));
    TA_rmse = TA_rmse + abs(min(cent2)-min([real_p2f1,real_p2f2]));
    TA_rmse = TA_rmse + abs(max(cent2)-max([real_p2f1,real_p2f2]));
end

xlabel('Estimated TA values (g)')
ylabel('Measured TA values (g)')
title('TA estimation performance for each leg across all interventions')
xlim([1 5])
ylim([1 5])
TA_rmse/(4*length(takes)) % final rmse

%% finding mean and std of measured TA to justify TA estimation rmse 6/15/22

takes = {'regular1', 'limp1', 'limp2', 'weight1', 'weight2', 'regular2'};
overall_mean = 0;
overall_std = 0;
for t = 1:length(takes)
    filename1 = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\',exp_subject,'\ProcessedData\ExcelData\both_', char(takes(t)),'_ta_p1_results.csv'];
    T1 = readtable(filename1);
    A1 = table2array(T1);
    real_TA1 = A1(2:end,2);
    overall_mean = overall_mean + mean(real_TA1);
    overall_std = overall_std + std(real_TA1);
    
    filename2 = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\',exp_subject,'\ProcessedData\ExcelData\both_', char(takes(t)),'_ta_p2_results.csv'];
    T2 = readtable(filename2);
    A2 = table2array(T2);
    real_TA2 = A2(2:end,2);
    overall_mean = overall_mean + mean(real_TA2);
    overall_std = overall_std + std(real_TA2);
end

avgmean = overall_mean/(2*length(takes)) % final mean
avgstd = overall_std/(2*length(takes)) % final mean

%% getting rmse on scaled GMM and original GMM results from estimated step times

takes = {'regular1', 'limp1', 'limp2', 'weight1', 'weight2', 'regular2'};
gmm_scaled_rmse = 0;
gmm_original_rmse = 0;
gmm_derived = 0;
for t = 1:length(takes)
    processedfilepath = ['C:\Users\katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\both_',char(takes(t)),'.mat'];
    load(processedfilepath)
    estimated_scaling_GMM_means = [0,0;0,0];
    estimated_original_GMM_means = [0,0;0,0];
    for p = 1:2 % for each person
        idx = find(clean_est_impacts(:,2) == p);
        estimated_impact_times = clean_est_impacts(idx,1)./Fs_pcb;

        % calculate all step times from impact times
        step_times = [];
        for x = 2:length(estimated_impact_times)
            curr_diff = estimated_impact_times(x) - estimated_impact_times(x-1);
            if curr_diff < 1.5 % not a turning point
                step_times(end+1) = curr_diff;
            end
        end
        GM = fitgmdist(transpose(step_times),2,'RegularizationValue',0.000001);
        proportion = GM.ComponentProportion;
        mu = GM.mu;
        [bigprop,~] = max(proportion);
        if (1-bigprop) < abs(0.5-bigprop)
            mean1 = mu(1) + (mu(2)-mu(1))*(1- (proportion(1))^2);
            mean2 = mu(2) + (mu(1)-mu(2))*(1- (proportion(2))^2);
        else
            mean1 = mu(1) + (mu(2)-mu(1))*abs(0.5-proportion(1));
            mean2 = mu(2) + (mu(1)-mu(2))*abs(0.5-proportion(2));
        end
        if p == 1
            estimated_scaling_GMM_means(1,:) = [min(mean1,mean2),max(mean1,mean2)];
            estimated_original_GMM_means(1,:) = [min(mu(1),mu(2)),max(mu(1),mu(2))];
        else
            estimated_scaling_GMM_means(2,:) = [min(mean1,mean2),max(mean1,mean2)];
            estimated_original_GMM_means(2,:) = [min(mu(1),mu(2)),max(mu(1),mu(2))];
        end
    end
    gmm_scaled_rmse = gmm_scaled_rmse + (estimated_scaling_GMM_means(1,1)-real_means(1,1))^2;
    gmm_scaled_rmse = gmm_scaled_rmse + (estimated_scaling_GMM_means(1,2)-real_means(1,2))^2;
    gmm_scaled_rmse = gmm_scaled_rmse + (estimated_scaling_GMM_means(2,1)-real_means(2,1))^2;
    gmm_scaled_rmse = gmm_scaled_rmse + (estimated_scaling_GMM_means(2,2)-real_means(2,2))^2;
    
    gmm_original_rmse = gmm_original_rmse + (estimated_original_GMM_means(1,1)-real_means(1,1))^2;
    gmm_original_rmse = gmm_original_rmse + (estimated_original_GMM_means(1,2)-real_means(1,2))^2;
    gmm_original_rmse = gmm_original_rmse + (estimated_original_GMM_means(2,1)-real_means(2,1))^2;
    gmm_original_rmse = gmm_original_rmse + (estimated_original_GMM_means(2,2)-real_means(2,2))^2;
    
    gmm_derived = gmm_derived + (estimated_scaled_means(1,1)-real_means(1,1))^2;
    gmm_derived = gmm_derived + (estimated_scaled_means(1,2)-real_means(1,2))^2;
    gmm_derived = gmm_derived + (estimated_scaled_means(2,1)-real_means(2,1))^2;
    gmm_derived = gmm_derived + (estimated_scaled_means(2,2)-real_means(2,2))^2;
end

gmm_scaled_rmse = sqrt(gmm_scaled_rmse/(4*length(takes)))
gmm_original_rmse = sqrt(gmm_original_rmse/(4*length(takes)))
gmm_derived = sqrt(gmm_derived/(4*length(takes)))


