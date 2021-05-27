%% load data 5/26/21
clear all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\';

datestr = '05-09-2021_15-08-29';
load([data_root_katie, 'DelsysFSRs\Run_number_110_', datestr])
load([data_root_katie, 'ImmersionFloorSensors\data_', datestr])
T = readtable([data_root_katie, 'MotionCapture\', datestr]);

% clip accel data from square triggers in fsr
fsrT = Time(1,:);
fsrD = Data;
[clipped_time, clipped_data] = clip_to_fsr(times, datas);

% convert and clip stomp mocap csv
[mocapT, mocapL, mocapR] = convertMocap(T);

%% todo****** extract footfall locations from these arrays

%% clip stomp fsr and pcb 

[maxv,~] = max(clipped_data(1:50000,2));
datalen = length(clipped_data);
buffer = 0.5;
for i = 1:datalen
   if clipped_data(i, 2) > maxv + buffer % always stomp next to S2
       break
   end
end

window = clipped_data(i-3000:i+10000,2); % approx. before and after impact
arrival_idx = aic_pick(window, 'to_peak')+i-3000;
clipped_time(arrival_idx) % this is the time pcb thinks impact happened




%% clip end of data

% [minv,~] = min(clipped_data(end-100:end,:));
[maxv,~] = max(clipped_data(end-100:end,:));

datalen = length(clipped_data);
for i = 1:datalen
    if any(clipped_data(datalen+1-i,:) > maxv)
        break
    end
end

buffer = 100;
lastindex = datalen+1-i + buffer;

clipped_data = clipped_data(0:lastindex,:);
clipped_time = clipped_time(0:lastindex);

lastindexfsr = findTindex(clipped_time(end),fsrT);
fsrT = fsrT(0:lastindexfsr);
fsrD = fsrD(:,0:lastindexfsr);

lastindexmc = findTindex(clipped_time(end),mocapT);
mocapL = mocapL(0:lastindexmc,:);
mocapR = mocapR(0:lastindexmc,:);
mocapT = mocapT(0:lastindexmc,:);

%% process mocap csv

mocapT = T(:,2);
T = T(:,3:end);

Lheel = 0;
Lmid = 1; 
Ltoe = 2;
Rheel = 3;
Rmid = 4;
Rtoe = 5;
S1 = 6;
S2 = 7;
S3 = 8;
S4 = 9;

S1loc = [T(1,S1*3 +1),T(1,S1*3 + 2),T(1,S1*3 + 3)];
S2loc = [T(1,S2*3 +1),T(1,S2*3 + 2),T(1,S2*3 + 3)];
S3loc = [T(1,S3*3 +1),T(1,S3*3 + 2),T(1,S3*3 + 3)];
S4loc = [T(1,S4*3 +1),T(1,S4*3 + 2),T(1,S4*3 + 3)];

%% sync FSR and pcb
invtrig = -1 .* datas(:,5);
edges = [0,0,0,0];
for n = 1:4
    [~,i] = max(invtrig);
    edges(n) = i;
    invtrig(i - 10000:i+10000) = -5;
end
edges = sort(edges);
startI = edges(1);
endI = edges(3);

accel = datas(startI:endI,1:4);
startT = times(startI);
endT = times(endI);
accelT = linspace(0,endT-startT,endI-startI+1);

% accel and accelT are synced to FSR (0 time are lined up)

%% process FSR
% double check later!
% FSR1: A outside ball, B toe C heel D inside ball
% FSR2: A toe, B heel, C inside ball, D outside ball

% figure;
% hold on
% for i = 1:4
%     subplot(4,1,i)
%     plot(Data(i+4,:))
% end

fsrT = Time(1,:);

Rballi  = Data(1,:);
Rtoe    = Data(2,:);
Rheel   = Data(3,:);
Rballo  = Data(4,:);

Ltoe    = Data(5,:);
Lheel   = Data(6,:);
Lballi  = Data(7,:);
Lballo  = Data(8,:);

figure
hold on
for i = 1:8
    subplot(8,1,i)
    plot(Data(i,10000:15000))
end


fsr = [Rheel;Rballi;Rballo;Rtoe;Lheel;Lballi;Lballo;Ltoe];

%% walking segments extract

wWidth = 1000; % number of datapoints in a typical impact window width
[~,seg1] = max(accel(:,2));
[~,seg2] = max(accel(:,3));
[~,seg3] = max(accel(seg1+wWidth:end,2));
seg3 = seg3 + seg1 + wWidth;

window1 = accel(seg1-wWidth:seg1+wWidth,2);
window2 = accel(seg2-wWidth:seg2+wWidth,3);
window3 = accel(seg3-wWidth:seg3+wWidth,2);

arrival1 = aic_pick(window1, 'to_peak') + seg1 - wWidth;
arrival2 = aic_pick(window2, 'to_peak') + seg2 - wWidth;
arrival3 = aic_pick(window3, 'to_peak') + seg3 - wWidth;

start1 = arrival1+wWidth;
end1 = arrival2-wWidth;
start2 = arrival2+wWidth;
end2 = arrival3-wWidth;

% final accelerometer segments
accelW1 = accel(start1:end1,:);
accelT1 = accelT(start1:end1);
accelW2 = accel(start2:end2,:);
accelT2 = accelT(start2:end2);

startT1 = accelT1(1);
endT1 = accelT1(end);
startT2 = accelT2(1);
endT2 = accelT2(end);

fsr_startI1 = (find(fsrT>startT1,1));
fsr_endI1 = (find(fsrT>endT1,1));
fsr_startI2 = (find(fsrT>startT2,1));
fsr_endI2 = (find(fsrT>endT2,1));

% final fsr segments
fsrT1 = fsrT(fsr_startI1:fsr_endI1);
fsrT2 = fsrT(fsr_startI2:fsr_endI2);
fsrW1 = fsr(:,fsr_startI1:fsr_endI1);
fsrW2 = fsr(:,fsr_startI2:fsr_endI2);

%% process pcb data

impactN = 1;
impactWidth = 10000; % number of indeces an impact lasts in data
mass = m; % mass of object if dropped

arrival_idx = zeros(n,sn);
peak_idx = zeros(n,sn);
peak_mag = zeros(n,sn);
% prominence = zeros(n,sn);
% echo = zeros(n,sn);
% area_under_curve = zeros(n,sn);

datas = filt_data;

for j = 1:sn
    for i = 1:impactN
        [~, maxidx] = max(datas(:,j));
        peak_idx(i,j) = maxidx;
        min_idx = maxidx - impactWidth*0.2;
        max_idx = min_idx + impactWidth;
        datas(min_idx:max_idx,j) = 0;
    end
end

%% making an example plot 5/12/21

% figure
% subplot(8,1,1)
% plot(fsrT1,fsrW1(1,:))
% subplot(8,1,2)
% plot(fsrT1,fsrW1(4,:))
% % subplot(8,1,3)
% % plot(accelT1,accelW1(:,1))
% % subplot(8,1,4)
% % plot(accelT1,accelW1(:,2))
% subplot(8,1,5)
% plot(accelT1,accelW1(:,3))
% subplot(8,1,6)
% plot(accelT1,accelW1(:,4))
% subplot(8,1,7)
% plot(fsrT1,fsrW1(5,:))
% subplot(8,1,8)
% plot(fsrT1,fsrW1(8,:))

% fsr = [Rheel;Rball1;Rball2;Rtoe;Lheel;Lball1;Lball2;Ltoe];

figure;
subplot(6,1,1)
plot(fsrT1(4234:4857),fsrW1(1,4234:4857))
hold on
xline(24.5237, 'b')
xline(25.056, 'r')
xline(25.2392, 'g')
xlabel('Time (s)')
title('Heel strike FSR')
subplot(6,1,2)
plot(fsrT1(4234:4857),fsrW1(4,4234:4857))
hold on
xline(24.5237, 'b')
xline(25.056, 'r')
xline(25.2392, 'g')
xlabel('Time (s)')
title('Toe strike FSR')
subplot(6,1,3)
plot(fsrT1(4234:4857),fsrW2(1,4234:4857))
hold on
xline(24.5237, 'b')
xline(25.056, 'r')
xline(25.2392, 'g')
xlabel('Time (s)')
title('Heel strike FSR')
subplot(6,1,4)
plot(fsrT1(4234:4857),fsrW2(4,4234:4857))
hold on
xline(24.5237, 'b')
xline(25.056, 'r')
xline(25.2392, 'g')
xlabel('Time (s)')
title('Toe strike FSR')
subplot(6,1,5)
plot(accelT1(104535:119895),accelW1(104535:119895,3))
hold on
xline(24.5237, 'b')
xline(25.056, 'r')
xline(25.2392, 'g')
xlabel('Time (s)')
title('Sensor 3 accelerometer')
subplot(6,1,6)
plot(accelT1(104535:119895),accelW1(104535:119895,4))
hold on
xline(24.5237, 'b')
xline(25.056, 'r')
xline(25.2392, 'g')
xlabel('Time (s)')
title('Sensor 4 accelerometer')

%% example plot 5/13/21

figure;
subplot(8,1,1)
plot(fsrT1(4234:4857),fsrW1(1,4234:4857))
hold on
xline(24.5237, 'b')
xline(25.056, 'r')
xline(25.2392, 'g')
xlabel('Time (s)')
title('Right Heel')
subplot(8,1,2)
plot(fsrT1(4234:4857),fsrW1(4,4234:4857))
hold on
xline(24.5237, 'b')
xline(25.056, 'r')
xline(25.2392, 'g')
xlabel('Time (s)')
title('Right Toe')
subplot(8,1,3)
plot(fsrT1(4234:4857),fsrW2(1,4234:4857))
hold on
xline(24.5237, 'b')
xline(25.056, 'r')
xline(25.2392, 'g')
xlabel('Time (s)')
title('Left Heel')
subplot(8,1,4)
plot(fsrT1(4234:4857),fsrW2(4,4234:4857))
hold on
xline(24.5237, 'b')
xline(25.056, 'r')
xline(25.2392, 'g')
xlabel('Time (s)')
title('Right Toe')

subplot(8,1,5)
plot(accelT1(104535:119895),accelW1(104535:119895,1))
hold on
xline(24.5237, 'b')
xline(25.056, 'r')
xline(25.2392, 'g')
xlabel('Time (s)')
title('Sensor 1 accelerometer')
subplot(8,1,6)
plot(accelT1(104535:119895),accelW1(104535:119895,2))
hold on
xline(24.5237, 'b')
xline(25.056, 'r')
xline(25.2392, 'g')
xlabel('Time (s)')
title('Sensor 2 accelerometer')
subplot(8,1,7)
plot(accelT1(104535:119895),accelW1(104535:119895,3))
hold on
xline(24.5237, 'b')
xline(25.056, 'r')
xline(25.2392, 'g')
xlabel('Time (s)')
title('Sensor 3 accelerometer')
subplot(8,1,8)
plot(accelT1(104535:119895),accelW1(104535:119895,4))
hold on
xline(24.5237, 'b')
xline(25.056, 'r')
xline(25.2392, 'g')
xlabel('Time (s)')
title('Sensor 4 accelerometer')
%%
peak_idx = sort(peak_idx, 'ascend');

for j = 1:sn
    for i = 1:n
        pk_i = peak_idx(i,j);
        min_idx = pk_i - impactWidth*0.2;
        max_idx = min_idx + impactWidth;
        window = filt_data(min_idx:max_idx,j);
        arrival_idx(i,j) = aic_pick(window, 'to_peak')+min_idx;
        peak_mag(i,j) = max(window);
        % prominence calc
        % echo calc
        % area under curve calc
    end
end

figure;
for i = 1:sn
    hold on
    subplot(sn,1,i)
    title(append('Filtered data at ', loc_names(i)))
    xlabel('Impact number')
    ylabel('Volts (V)')
    hold on
    plot(filt_data(:,i))
    plot(arrival_idx(:,i),0,'bx')
    plot(peak_idx(:,i),peak_mag(:,i),'ro')
end

% save([data_root_katie, filename], 'raw_data','loc_names','Fs','filt_data','impactN','arrival_idx','peak_idx','peak_mag','mass', '-append')
   
%%
fsr = Data(1,:);
Time = Time(1,:);

for i = 1:length(fsr)
    if fsr(i) < 0.02
        fsr(i) = 0;
    end
end




