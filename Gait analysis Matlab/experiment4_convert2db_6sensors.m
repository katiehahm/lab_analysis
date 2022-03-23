%% 3/11/22 experiment 4 with Praneeth processing
clear all
close all
filepath = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Praneeth experiment 3_11_22\both_';
% ################# change ########################
intervention = 'towards1'; % regular1 brace1 brace2 weight1 weight2 regular2
% #################################################

% load 3 datasets
load([filepath, intervention])
load([filepath, intervention, '_fsr'])
T = readtable([filepath, intervention, '_mocap']);

% constants
sensorN = 6;
Fs_pcb = 12800;
Fs_fsr = Fs(1);
Fs_acc = Fs(2);
TriggerRow = 9; % row of Data that contains trigger change #####################
MCkeySet = {'Lx1','Ly1','Lz1','Rx1','Ry1','Rz1','Lx2','Ly2','Lz2','Rx2','Ry2','Rz2'};
MCvalueSet = [3 4 5 6 7 8 9 10 11 12 13 14]; % order of Mocap left/right based on keyset
Mmocap = containers.Map(MCkeySet, MCvalueSet);

[mocapT, mocapL1, mocapR1, mocapL2, mocapR2] = convertMocap2pp(T, Mmocap);
allmocap = cat(3, mocapL1,mocapR1,mocapL2,mocapR2);

[accData, accTime, fsrData, fsrTime] = clip_accelfsr_fromMocap2pp(Data, Time, Fs, TriggerRow);
[pcbData, pcbTime] = clip_pcb_fromMocap(datas, times, sensorN);

% check this is correct. These numbers should be nearly the same:
mocapT(end)
fsrTime(end)
pcbTime(end)
accTime(end)

%% delete end of data when walking off
figure; plot(mocapT, mocapR1(:,1),'b.')
figure; plot(pcbTime, pcbData)

%% cont.
last_time = 87.13;
lastidx_pcb = findTindex(last_time,pcbTime);
lastidx_fsr = findTindex(last_time,fsrTime);
lastidx_acc = findTindex(last_time,accTime);
lastidx_mocap = findTindex(last_time,mocapT);

fsrData = fsrData(1:lastidx_fsr,:);
fsrTime = fsrTime(1:lastidx_fsr);
pcbTime = pcbTime(1:lastidx_pcb);
pcbData = pcbData(1:lastidx_pcb,:);
mocapT = mocapT(1:lastidx_mocap);
allmocap = allmocap(1:lastidx_mocap,:,:);
accTime = accTime(1:lastidx_acc);
accData = accData(1:lastidx_acc,:);

mocapT(end)
fsrTime(end)
pcbTime(end)
accTime(end)

%%
% clean and filter pcb data
filt_pcbD = lpf_data(pcbData);

% finding footfalls based on fsr heel data
L_dist = 180; % min distance between heel peaks
R_dist = 180;
min_threshL = 25; % value between peaks threshold, if not lower than this then omit peak
min_threshR = 25;

% min distance between heel peaks
dist = [180, 180, 180, 180]; % L1, R1, L2, R2
min_thresh = [25, 25, 25, 25]; % L1, R1, L2, R2

impacts = findimpacts_fsr_accel2pp(fsrTime, mocapT, allmocap, fsrData,dist,min_thresh);

%% fix small errors in impacts 11/5/21
heel_start_wrong = [2846,7679,7948]; % these need to be same length
heel_start_right = [2834,7662,7935];

heel_pk_wrong = []; % index, these need to be same length
heel_pk_right = [];

delete_pks = [18895]; % index of peaks to delete

impacts = manual_fix_fsr2pp(impacts,fsrData,heel_start_wrong,heel_start_right,heel_pk_wrong,heel_pk_right,delete_pks);

%% pcb adjust params to capture whole impact 11/5/21
% look at where the impacts are and adjust window_width and offset
% parameters in the next section

close all
secs = 20; % plot 10 secs of data
figure;
for i = 1:6
    subplot(6,1,i)
    plot(pcbTime(1:secs*Fs_pcb), filt_pcbD(1:secs*Fs_pcb,i))
    hold on
    indeces = find(fsrTime(impacts(:,1)) < secs);
    plot(fsrTime(impacts(indeces,1)),0,'r.','MarkerSize',10)
end

%% pcb and mocap extract

window_width = Fs_pcb*0.3; % shaking lasts < N seconds
offset = 0.1; % start window N behind fsr start
[arrival_idx, peak_idx, peak_mag] = findimpacts_pcb(window_width,offset,impacts,fsrTime,pcbTime,filt_pcbD,sensorN,true);

%% mocap extract
[coordinates] = findimpacts_mocap2pp(impacts,fsrTime,mocapT,allmocap,true);

%% extract accel data
% gets the peak of abs value of accelerometer values for x,y,z directions
% 4th col in variable 'acc_pks' indicates whether it is LorR leg.
window = 0.3; % seconds
offset = 0.25;

[acc_pks,acc_pk_idx] = find_accel_impacts2pp(impacts,fsrTime,window,offset,accTime,accData,fsrData);

%% saving data to matlab
processedfilepath = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Praneeth experiment 3_11_22\ProcessedData\';
filename = [processedfilepath, 'both_',intervention];

save(filename,'pcbTime','filt_pcbD','arrival_idx','peak_idx','peak_mag','Data','datas','pcbData', ...
    'fsrTime','fsrData','impacts','Fs_pcb','Fs_fsr','Fs_acc', ...
    'acc_pks','acc_pk_idx','accTime','accData',...
    'mocapT','allmocap','coordinates','sensorN')
disp(append("Saved as ", filename))

%% make plot just for com mtg 2/9/22

% for i = 1:4
%     figure;
%     p1 = plot(pcbTime, filt_pcbD(:,i),'k');
%     hold on
%     person1_idx = find(coordinates(:,5) == 1 | coordinates(:,5) == 2);
%     p2 = plot(mocapT(coordinates(person1_idx,4))+0.08,0,'r.','MarkerSize',20);
% %     for n = 1:length(person1_idx)
% %         xline(mocapT(coordinates(person1_idx(n),4))+0.08,'--r')
% %     end
%     person2_idx = find(coordinates(:,5) == 3 | coordinates(:,5) == 4);
%     p3 = plot(mocapT(coordinates(person2_idx,4))+0.08,0,'g.','MarkerSize',20);
% %     for n = 1:length(person2_idx)
% %         xline(mocapT(coordinates(person2_idx(n),4))+0.08,'--b')
% %     end
%     h = [p1(1);p2(1);p3(1)];
%     legend(h, 'Accelerometer','Person 1','Person 2')
%     xlabel('Time (s)')
%     ylabel('Volts (V)')
% end
    
figure;
p1 = plot(pcbTime(441602:490241), filt_pcbD(441602:490241,1),'k');
hold on
person1_idx = find(coordinates(:,5) == 1 | coordinates(:,5) == 2);
reds = mocapT(coordinates(person1_idx,4));
reds_idx = find(reds > 34.5);
reds_idx2 = find(reds > 38.3);
p2 = plot(reds(reds_idx(1):reds_idx2(1)-1)+0.08,0,'r.','MarkerSize',20);
person2_idx = find(coordinates(:,5) == 3 | coordinates(:,5) == 4);
greens = mocapT(coordinates(person2_idx,4));
greens_idx = find(greens > 34.5);
greens_idx2 = find(greens > 38.3);
p3 = plot(greens(greens_idx(1):greens_idx2(1)-1)+0.08,0,'g.','MarkerSize',20);
h = [p1(1);p2(1);p3(1)];
legend(h, 'Accelerometer','Person 1','Person 2')
xlabel('Time (s)')
ylabel('Volts (V)')
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14)






