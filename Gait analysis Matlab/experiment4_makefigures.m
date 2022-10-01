%% single person data 7/13/22
load('C:\Users\katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Alex 4\alone_regular.mat')
figure; plot(times,datas(:,1))
xlim([54 60])
ylim([-0.015 0.015])
xlabel('Time (s)')
ylabel('Volts (V)')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)


%% change data in a figure without source code 7/8/22
fig = gcf;
a = get(gca,'Children');
ydata = get(a,'YData');
xdata = get(a,'XData');

X = cell2mat(xdata);
Y = cell2mat(ydata);
X = X(3:end);
Y = Y(3:end);

figure;
meas = find(Y == 0);
detect = find(Y == 0.5);
plot(X(meas),Y(meas),'bx','MarkerSize',12)
hold on
plot(X(detect),Y(detect),'rx','MarkerSize',12)
legend('Measured impacts','Detected impacts')
ylim([-0.5 0.9])
xlim([47.3 52.2])
title('Detected impacts from all sensors')
xlabel('Time (s)')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)

%% exponent pres, raw signal vs scalogram vs sum of scalogram mag

load('C:\Users\katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Jenny 1\ProcessedData\both_regular1.mat')
starttime = 1.25;
lasttime = 3.25;

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

%% exponent pres, edit to get wiener signal

load('C:\Users\katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Jenny 1\ProcessedData\both_regular1.mat')
starttime = 1.25;
lasttime = 3.25;

% raw signal
figure;
pcbstarttime = starttime*Fs_pcb;
pcblasttime = lasttime*Fs_pcb;
pcbclip = wien_pcbD(pcbstarttime:pcblasttime,1);
plot(pcbTime(pcbstarttime:pcblasttime)-starttime,pcbclip)
xlabel('Time (s)')
ylabel('Volts (V)')
title('Wiener filter accelerometer signal')
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
b = bar(X,GMMresults);
legend('Measured mean subj 1','Estimated mean subj 1','Measured mean subj 2','Estimated mean subj 2')
ylabel('Step time (s)')
ylim([0.3,0.75])
title('Estimated and measured means for one experiment')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)
b(1).FaceColor = [1 0 0];
b(2).FaceColor = [1 0.5 0.5];
b(3).FaceColor = [0 0 1];
b(4).FaceColor = [0.5 0.5 1];

%% GMM results all subjects

takes = {'regular1', 'limp1', 'limp2', 'weight1', 'weight2', 'regular2'};
subj = {'Jenny 1','Uriel 2','April 3'};
% subj = {'Uriel 2'};

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
GMM(1,:) = [];
figure;
one2one = linspace(-4,4,1000);
plot(one2one,one2one,'k--','LineWidth',0.5)
hold on
plot(GMM(:,1),GMM(:,2),'b.','MarkerSize',13)
xlabel('Estimated step time (s)')
ylabel('Measured step time (s)')
title('Estimated and measured step times for all subjects')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)
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
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)



%% localization results matlab fig

filename = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Jenny 1\ProcessedData\ExcelData\alltakes_localization_p1_results.csv'];
T = readtable(filename);
A = table2array(T);
estlocs = A(:,2);
reallocs = A(:,1);
figure;
one2one = linspace(-4,4,1000);
plot(one2one,one2one,'k--','LineWidth',2)
hold on
plot(estlocs, reallocs,'b.')
xlabel('Estimated locations (m)')
ylabel('Measured locations (m)')
title('Localization performance for one example subject')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)

%% GMM good and bad plot

r1 = normrnd(0.5,0.1,[1,500]);
r2 = normrnd(0.52,0.1,[1,500]);

figure;
histogram(r1,'FaceAlpha',0.5,'EdgeAlpha',0.5,'FaceColor','#0072BD')
hold on
% hold on
% y = 0.3:0.1:0.7;
% mu = 0.5;
% sigma = 0.1;
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'LineWidth',1.5)
histogram(r2,'FaceAlpha',0.5,'EdgeAlpha',0.5,'FaceColor','#D95319')
% hold on
% y = 0.3:0.1:0.7;
% mu = 0.52;
% sigma = 0.1;
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'LineWidth',1.5)
xline(0.5,'LineWidth',2,'Color','b')
xline(0.52,'LineWidth',2,'Color','r')
legend('Leg 1 ST','Leg 2 ST','Leg 1 avg ST','Leg 2 avg ST')
title('Step time (ST) distribution derived GMM results')
ylabel('Occurances')
xlabel('Step time (s)')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)

r_all = [r1,r2];
GMModel = fitgmdist(r_all.',2); % this is commented bc got a good split 
f1 = normrnd(GMModel.mu(1),0.1,[1,round(1000*GMModel.ComponentProportion(1))]);
f2 = normrnd(GMModel.mu(2),0.1,[1,round(1000*GMModel.ComponentProportion(2))]);
figure;
histogram(f1,'FaceAlpha',0.5,'EdgeAlpha',0.5,'FaceColor','#0072BD')
hold on
% y = 0.3:0.1:0.7;
% mu = GMModel.mu(1);
% sigma = 0.1;
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'LineWidth',1.5)
histogram(f2,'FaceAlpha',0.5,'EdgeAlpha',0.5,'FaceColor','#D95319')
% y = 0.3:0.1:0.7;
% mu = GMModel.mu(2);
% sigma = 0.1;
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'LineWidth',1.5)
xline(GMModel.mu(1),'LineWidth',2,'Color','b')
xline(GMModel.mu(2),'LineWidth',2,'Color','r')

legend('Leg 1 ST','Leg 2 ST','Leg 1 avg ST','Leg 2 avg ST')
title('Step time (ST) distribution original GMM results')
ylabel('Occurances')
xlabel('Step time (s)')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)

%% TA estimation for one exp (copied from _processing5.m)
takes = {'regular1', 'limp1', 'limp2', 'weight1', 'weight2', 'regular2'};
exp_subject = 'Jenny 1';
TA_rmse = 0;
figure;
one2one = linspace(1,5,1000);
plot(one2one,one2one,'k--','LineWidth',0.5)
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
    foot_labels = correct_coords(:,3);
    p1f1_idx = find(foot_labels == 1);
    p1f2_idx = find(foot_labels == 2);
    p2f1_idx = find(foot_labels == 3);
    p2f2_idx = find(foot_labels == 4);
    real_p1f1 = mean(correct_ta(p1f1_idx));
    real_p1f2 = mean(correct_ta(p1f2_idx));
    real_p2f1 = mean(correct_ta(p2f1_idx));
    real_p2f2 = mean(correct_ta(p2f2_idx));
    
    plot(min(cent1),min([real_p1f1,real_p1f2]),'b.','MarkerSize',15)
    plot(max(cent1),max([real_p1f1,real_p1f2]),'b.','MarkerSize',15)
    plot(min(cent2),min([real_p2f1,real_p2f2]),'b.','MarkerSize',15)
    plot(max(cent2),max([real_p2f1,real_p2f2]),'b.','MarkerSize',15)
    
    % calculate rmse
    TA_rmse = TA_rmse + abs(min(cent1)-min([real_p1f1,real_p1f2]));
    TA_rmse = TA_rmse + abs(max(cent1)-max([real_p1f1,real_p1f2]));
    TA_rmse = TA_rmse + abs(min(cent2)-min([real_p2f1,real_p2f2]));
    TA_rmse = TA_rmse + abs(max(cent2)-max([real_p2f1,real_p2f2]));
end

xlabel('Estimated TA values (g)')
ylabel('Measured TA values (g)')
title('TA estimation performance for one experiment')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)
xlim([1 5])
ylim([1 5])
TA_rmse/(4*length(takes)) % final rmse




