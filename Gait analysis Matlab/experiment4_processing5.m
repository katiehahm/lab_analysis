% makes overall csv for alltakes localization, TA estimation, kmeans
% 6/10/22
%% making a localization csv for all takes

takes = {'regular1', 'limp1', 'limp2', 'weight1', 'weight2', 'regular2'};
newA = zeros(1,53);
person = '2';
for t = 1:length(takes)
    filename = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\ExcelData\both_', char(takes(t)),'_localization_p',person','_withta.csv'];
    T = readtable(filename);
    A = table2array(T);
    newA = [newA;A];
end
newA(1,:) = []; % initialization

filename = ['C:/Users/Katie/Dropbox (MIT)/Lab/Analysis/Experiment4/April 3/ProcessedData/ExcelData/alltakes_localization_p',person','_withta.csv'];
writematrix(newA,filename)

%% make TA csv with overlapping impacts 5/26/22 using alltakes localization results

filename1 = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\ExcelData\alltakes_localization_p1_results.csv';
T = readtable(filename1);
loc_results1 = table2array(T);
filename2 = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\ExcelData\alltakes_localization_p2_results.csv';
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
        filename = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\ExcelData\both_', char(takes(t)),'_localization_p',person','_withta.csv'];
        T = readtable(filename);
        A = table2array(T);
        numrows = length(A(:,1));
        if p == 1
            pred_locs = loc_results1(:,2); % second col has est values
        else
            pred_locs = loc_results2(:,2);
        end
        locs = pred_locs(rowcount:rowcount + numrows - 1);
        featureV = zeros(1,31);
        for i = 1:numrows
            xcoord = locs(i);
            dist1 = sqrt( (xcoord-s1(1)).^2 );
            dist2 = sqrt( (xcoord-s2(1)).^2 );
            dist3 = sqrt( (xcoord-s3(1)).^2 );
            dist4 = sqrt( (xcoord-s4(1)).^2 );
            dist5 = sqrt( (xcoord-s5(1)).^2 );
            dist6 = sqrt( (xcoord-s6(1)).^2 );
            feature = [A(i,4),dist1,dist2,dist3,dist4,dist5,dist6,A(i,5:28)];
            featureV(end+1,:) = feature;
        end
        featureV(1,:) = []; % initialization
        newfilename = ['C:/Users/Katie/Dropbox (MIT)/Lab/Analysis/Experiment4/April 3/ProcessedData/ExcelData/both_', char(takes(t)),'_ta_p',person,'.csv'];
        writematrix(featureV,newfilename)
        rowcount = rowcount + numrows;
    end
end

% then run exp4_TAestimation.py

%% doing kmeans on TA estimation results all takes 5/24/22

takes = {'regular1', 'limp1', 'limp2', 'weight1', 'weight2', 'regular2'};
newA = zeros(1,53);
TA_rmse = 0;
figure;
hold on
for t = 1:length(takes)
    filename1 = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\ExcelData\both_', char(takes(t)),'_ta_p1_results.csv'];
    T1 = readtable(filename1);
    A1 = table2array(T1);
    est_TA1 = A1(2:end,3);
    [~,cent1] = kmeans(est_TA1,2);
    
    filename2 = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\ExcelData\both_', char(takes(t)),'_ta_p2_results.csv'];
    T2 = readtable(filename2);
    A2 = table2array(T2);
    est_TA2 = A2(2:end,3);
    [~,cent2] = kmeans(est_TA2,2);
    
    load(['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\April 3\ProcessedData\both_',char(takes(t))])
    foot_labels = correct_coords(:,3);
    p1f1_idx = find(foot_labels == 1);
    p1f2_idx = find(foot_labels == 2);
    p2f1_idx = find(foot_labels == 3);
    p2f2_idx = find(foot_labels == 4);
    real_p1f1 = mean(correct_ta(p1f1_idx));
    real_p1f2 = mean(correct_ta(p1f2_idx));
    real_p2f1 = mean(correct_ta(p2f1_idx));
    real_p2f2 = mean(correct_ta(p2f2_idx));
    
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