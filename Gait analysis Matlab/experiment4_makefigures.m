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