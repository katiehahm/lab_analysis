% exponent pres, raw signal vs scalogram vs sum of scalogram mag

load('C:\Users\katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Jenny 1\ProcessedData\both_regular1.mat')
starttime = 1;
lasttime = 3;

% raw signal
figure;
pcbstarttime = starttime*Fs_pcb;
pcblasttime = lasttime*Fs_pcb;
pcbclip = pcbData(pcbstarttime:pcblasttime,1);
plot(pcbTime(pcbstarttime:pcblasttime)-starttime,pcbclip)
xlabel('Time (s)')
ylabel('Volts (V)')
title('Raw accelerometer signal')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)
xlim([0 lasttime-starttime])

% scalogram
figure;
cwt(pcbclip,Fs_pcb)
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)
[wt,f] = cwt(pcbclip,Fs_pcb); % uses default Morse wavelet
valid_f_idx = find(f < 350 & f > 100);
cwt_freq = f(valid_f_idx);
cwt_mag = abs(wt(valid_f_idx,:));
sum_cwt = sum(cwt_mag,1);
sum_smooth_cwt = movmean(sum_cwt, 800);

%% sum scalogram
figure;
cwt_time = linspace(starttime,lasttime,length(sum_smooth_cwt));
plot(cwt_time-starttime,sum_smooth_cwt)
hold on
fsrstarttime = starttime*Fs_fsr;
fsrlasttime = lasttime*Fs_fsr;
impactidx = find(impacts(:,1) >= fsrstarttime & impacts(:,1) <= fsrlasttime);
impacttimes = impacts(impactidx,1)/Fs_fsr;
plot(impacttimes-starttime,0,'bx','MarkerSize',12)
hold on
idx = find(final_estimates(:,1) >= pcbstarttime & final_estimates(:,1) <= pcblasttime);
detect_impact = final_estimates(idx,1)./Fs_pcb;
plot(detect_impact-starttime,0,'rx','MarkerSize',12)
xlabel('Time (s)')
ylabel('Magnitudes')
title('Sum of CWT scalogram magnitudes')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)
xlim([0 lasttime-starttime])
legend('CWT magnitude sum')
plot([NaN NaN], [NaN NaN], 'x','Color', "blue", 'DisplayName', "Measured impact times")
plot([NaN NaN], [NaN NaN], 'x','Color', "red", 'DisplayName', "Detected impact times")

%% impacts detection, DFS sequence, recusive edit 5/3/22

starttime = 18.6;
lasttime = 21.2;
pcbstarttime = starttime*Fs_pcb;
pcblasttime = lasttime*Fs_pcb;
fsrstarttime = starttime*Fs_fsr;
fsrlasttime = lasttime*Fs_fsr;
impactidx = find(impacts(:,1) >= fsrstarttime & impacts(:,1) <= fsrlasttime);
impacttimes = impacts(impactidx,1)/Fs_fsr;
idx = find(final_estimates(:,1) >= pcbstarttime & final_estimates(:,1) <= pcblasttime);
detect_impact = final_estimates(idx,1)./Fs_pcb;
figure;
plot(detect_impact,1,'r*')
hold on
plot(impacttimes,0,'b*')
ylim([-1 2])
xlabel('Time (s)')
% legend('Detected impact times','Measured impact times')
title('Detected impacts compared with ground truth')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)
legend('Detected')
plot([NaN NaN], [NaN NaN], '*','Color', "blue", 'DisplayName', "Measured")

% manually make plot with DFS sequence
detect_impact = sort(detect_impact);
dfs_x = [detect_impact(1),detect_impact(3),detect_impact(5),detect_impact(7)];
dfs_o = [detect_impact(2),detect_impact(4),detect_impact(6),detect_impact(8)];
figure;
plot(dfs_x,1,'rx')
hold on
plot(dfs_o,1,'ro')
curr_ID = ID_labels(impactidx);
plot(impacttimes(curr_ID == 1),0,'bo')
plot(impacttimes(curr_ID == 2),0,'bx')
ylim([-1 2.5])
xlim([18.5 21.5])
xlabel('Time (s)')
% legend('Detected impact times','Measured impact times')
title('Initial person labels')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)
legend('Estimated person 1')
plot([NaN NaN], [NaN NaN], 'o','Color', "red", 'DisplayName', "Estimated person 2")
plot([NaN NaN], [NaN NaN], 'x','Color', "blue", 'DisplayName', "Measured person 1")
plot([NaN NaN], [NaN NaN], 'o','Color', "blue", 'DisplayName', "Measured person 2")

%% recursive sequence 1
dfs_x = [detect_impact(1),detect_impact(3),detect_impact(4),detect_impact(6),detect_impact(8)];
dfs_o = [detect_impact(2),detect_impact(3),detect_impact(5),detect_impact(7)];
figure;
plot(dfs_x,1,'rx')
hold on
plot(dfs_o,1,'ro')
curr_ID = ID_labels(impactidx);
plot(impacttimes(curr_ID == 1),0,'bo')
plot(impacttimes(curr_ID == 2),0,'bx')
ylim([-1 2.5])
xlim([18.5 21.5])
xlabel('Time (s)')
% legend('Detected impact times','Measured impact times')
title('Recursive sequence 1')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)
legend('Estimated person 1')
plot([NaN NaN], [NaN NaN], 'o','Color', "red", 'DisplayName', "Estimated person 2")
plot([NaN NaN], [NaN NaN], 'x','Color', "blue", 'DisplayName', "Measured person 1")
plot([NaN NaN], [NaN NaN], 'o','Color', "blue", 'DisplayName', "Measured person 2")

%% recursive sequence 2
dfs_x = dfs_x(1:end-1);

figure;
plot(dfs_x,1,'rx')
hold on
plot(dfs_o,1,'ro')
curr_ID = ID_labels(impactidx);
plot(impacttimes(curr_ID == 1),0,'bo')
plot(impacttimes(curr_ID == 2),0,'bx')
ylim([-1 2.5])
xlim([18.5 21.5])
xlabel('Time (s)')
% legend('Detected impact times','Measured impact times')
title('Recursive sequence 2')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)
legend('Estimated person 1')
plot([NaN NaN], [NaN NaN], 'o','Color', "red", 'DisplayName', "Estimated person 2")
plot([NaN NaN], [NaN NaN], 'x','Color', "blue", 'DisplayName', "Measured person 1")
plot([NaN NaN], [NaN NaN], 'o','Color', "blue", 'DisplayName', "Measured person 2")

%% impacts detection, DFS sequence, recusive edit
% same as above but NO GROUND TRUTH 5/4/22

starttime = 18.6;
lasttime = 21.2;
pcbstarttime = starttime*Fs_pcb;
pcblasttime = lasttime*Fs_pcb;
fsrstarttime = starttime*Fs_fsr;
fsrlasttime = lasttime*Fs_fsr;
impactidx = find(impacts(:,1) >= fsrstarttime & impacts(:,1) <= fsrlasttime);
impacttimes = impacts(impactidx,1)/Fs_fsr;
idx = find(final_estimates(:,1) >= pcbstarttime & final_estimates(:,1) <= pcblasttime);
detect_impact = final_estimates(idx,1)./Fs_pcb;
figure;
plot(detect_impact,0,'r*')
hold on
ylim([-1 2])
xlabel('Time (s)')
% legend('Detected impact times','Measured impact times')
title('Detected impact times')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)

% manually make plot with DFS sequence
detect_impact = sort(detect_impact);
dfs_x = [detect_impact(1),detect_impact(3),detect_impact(5),detect_impact(7)];
dfs_o = [detect_impact(2),detect_impact(4),detect_impact(6),detect_impact(8)];
figure;
plot(dfs_x,0,'rx')
hold on
plot(dfs_o,0,'ro')
ylim([-1 2])
xlim([18.5 21.5])
xlabel('Time (s)')
title('Initial person labels')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)
legend('Estimated person 1')
plot([NaN NaN], [NaN NaN], 'o','Color', "red", 'DisplayName', "Estimated person 2")

% recursive sequence 1
dfs_x = [detect_impact(1),detect_impact(3),detect_impact(4),detect_impact(6),detect_impact(8)];
dfs_o = [detect_impact(2),detect_impact(3),detect_impact(5),detect_impact(7)];
figure;
plot(dfs_x,0,'rx')
hold on
plot(dfs_o,0,'ro')
ylim([-1 2])
xlim([18.5 21.5])
xlabel('Time (s)')
title('Recursive sequence 1')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)
legend('Estimated person 1')
plot([NaN NaN], [NaN NaN], 'o','Color', "red", 'DisplayName', "Estimated person 2")

% recursive sequence 2
dfs_x = dfs_x(1:end-1);

figure;
plot(dfs_x,0,'rx')
hold on
plot(dfs_o,0,'ro')
ylim([-1 2])
xlim([18.5 21.5])
xlabel('Time (s)')
% legend('Detected impact times','Measured impact times')
title('Recursive sequence 2')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)
legend('Estimated person 1')
plot([NaN NaN], [NaN NaN], 'o','Color', "red", 'DisplayName', "Estimated person 2")

%% GMM results one example subject

takes = {'regular1', 'limp1', 'limp2', 'weight1', 'weight2', 'regular2'};
GMMresults = zeros(1,4);
for t = 1:length(takes)
    filename = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Uriel 2\ProcessedData\both_', char(takes(t))];
    load(filename)
    
    gmmarray = [real_means(1,1),estimated_scaled_means(1,1),real_means(1,2),estimated_scaled_means(1,2)];
    GMMresults(end+1,:) = gmmarray;
end

GMMresults(1,:) = []; % initialization
X = categorical({'Regular 1','Limp 1','Limp 2','Weight 1','Weight 2','Regular 2'});
X = reordercats(X,{'Regular 1','Limp 1','Limp 2','Weight 1','Weight 2','Regular 2'});

figure;
bar(X,GMMresults)
legend('Real mean 1','Estimated mean 1','Real mean 2','Estimated mean 2')
ylabel('Step time (s)')
title('Estimated and measured means of each leg for one subject')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 14)

%% GMM results all subjects

takes = {'regular1', 'limp1', 'limp2', 'weight1', 'weight2', 'regular2'};
subj = {'Jenny 1','Uriel 2'};

GMM = zeros(1,2); % [est, real]

for s = 1:length(subj)
    for t = 1:length(takes)
        filename = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\',char(subj(s)),'\ProcessedData\both_', char(takes(t))];
        load(filename)
        
        GMM(end+1,:) = [estimated_scaled_means(1,1),real_means(1,1)];
        GMM(end+1,:) = [estimated_scaled_means(1,2),real_means(1,2)];
        GMM(end+1,:) = [estimated_scaled_means(2,1),real_means(2,1)];
        GMM(end+1,:) = [estimated_scaled_means(2,2),real_means(2,2)];
    end
end

figure; 
plot(GMM(:,1),GMM(:,2),'b.')
xlabel('Estimated step time (s)')
ylabel('Measured step time (s)')
title('Estimated and measured means across all subjects')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 14)
xlim([0.35 0.9])
ylim([0.35 0.9])

sqrt(mean((GMM(:,1) - GMM(:,2)).^2))


%% variation in TA data 

Take = [];
Foot = [];
TA = [];
takes = {'regular1', 'limp1', 'limp2', 'weight1', 'weight2', 'regular2'};
for t = 1:length(takes)
    filename = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Jenny 1\ProcessedData\both_', char(takes(t))];
    load(filename)
    foot1 = find(impacts(:,4) == 3);
    foot2 = find(impacts(:,4) == 4);
    Take = [Take, ones(1,length(foot1) + length(foot2)).*t];
    Foot = [Foot, zeros(1,length(foot1)), ones(1,length(foot2))];
    TA = [TA, acc_pks(foot1,2).', acc_pks(foot2,2).'];
end

tbl = table(Foot, Take, TA);

figure;
boxchart(tbl.Take,tbl.TA,'GroupByColor',tbl.Foot,'JitterOutliers','off')
ylabel('Acceleration (g)')
legend('Right Foot','Left Root')
xticks([1,2,3,4,5,6]);
xticklabels({'regular1','brace1','brace2','weight1','weight2','regular2'})
title('Observed TA of each foot for one subject')
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14)


%% localization results matlab fig

filename = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Uriel 2\ProcessedData\ExcelData\alltakes_localization_p2_results.csv'];
T = readtable(filename);
A = table2array(T);
estlocs = A(:,2);
reallocs = A(:,1);
figure;
plot(estlocs, reallocs,'b.')
xlabel('Estimated locations (m)')
ylabel('Measured locations (m)')
title('Localization performance for one example subject')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 14)
