% load data 5/26/21
clear all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\';

datestr = '05-09-2021_15-08-29';
load([data_root_katie, 'DelsysFSRs\Run_number_110_', datestr])
load([data_root_katie, 'ImmersionFloorSensors\data_', datestr])
T = readtable([data_root_katie, 'MotionCapture\', datestr]);

% clip accel data from square triggers in fsr
Time = Time(1,:);
[clipped_time, clipped_data] = clip_to_fsr(times, datas);

% convert and clip stomp mocap csv
[mocapT, mocapL, mocapR] = convertMocap(T);

% todo****** extract footfall locations from these arrays

% clip stomp fsr and pcb 
[clip_stomp_data,clip_stomp_time,clip_fsrD,clip_fsrT] = clip_stomp(clipped_time,clipped_data,Time,Data);


% clip end of data
[pcbD,pcbT,fsrT,fsrD,mocapT,mocapR,mocapL] = clip_end(clip_stomp_data,clip_stomp_time,clip_fsrT,clip_fsrD,mocapT,mocapR,mocapL);

% clip 1 sec from start of data to filter out 
[pcbD,pcbT,fsrT,fsrD,mocapT,mocapR,mocapL] = clip_onesec(pcbD,pcbT,fsrT,fsrD,mocapT,mocapR,mocapL);

% map of fsr sensor location to array number
keySet = {'Rballi','Rtoe','Rheel','Rballo','Ltoe','Lheel','Lballi','Lballo'};
valueSet = [1 2 3 4 5 6 7 8];
Mfsr = containers.Map(keySet, valueSet);
% overall plot for visual check
plot_3data(pcbD,pcbT,fsrT,fsrD,mocapT,mocapR,mocapL,Mfsr)

% clean and filter pcb data
Fs = 12800;
filt_pcbD = lpf_data(pcbD);
filt_pcbD = noise_to_zero(filt_pcbD,Fs);

% Final workable data arrays:
% fsrD, fsrT, filt_pcbD, pcbT, mocapT, mocapR, mocapL


%% finding footfalls based on fsr heel data 6/7/21

figure;
plot(fsrD(:,Mfsr('Rheel')))
findpeaks(fsrD(:,Mfsr('Rheel')),'MinPeakProminence',10,'Annotate','extents','MinPeakDistance',0.5)
[pksR,locsR,widthsR,~] = findpeaks(fsrD(:,Mfsr('Rheel')),'MinPeakProminence',10,'Annotate','extents','MinPeakDistance',0.5); % 4th element is prominence

width_thresh = 400;
for i = 1:length(widthsR) % filter out wide peaks
    if widthsR(i) > width_thresh
        pksR(i) = 0;
        locsR(i) = 0;
    end
end
pksR = nonzeros(pksR);
locsR = nonzeros(locsR);

figure;
plot(fsrD(:,Mfsr('Lheel')))
findpeaks(fsrD(:,Mfsr('Lheel')),'MinPeakProminence',10,'Annotate','extents','MinPeakDistance',0.5)
[pksL,locsL,widthsL,~] = findpeaks(fsrD(:,Mfsr('Lheel')),'MinPeakProminence',10,'Annotate','extents','MinPeakDistance',0.5);

for i = 1:length(widthsL) % filter out wide peaks
    if widthsL(i) > width_thresh
        pksL(i) = 0;
        locsL(i) = 0;
    end
end
pksL = nonzeros(pksL);
locsL = nonzeros(locsL);

impacts = [locsR(1),pksR(1);
        locsL(1),pksL(1)];
for i = 2:length(pksR)
    btw_pk = fsrD(locsR(i-1):locsR(i),Mfsr('Rheel'));
    if any(btw_pk < 1)
        impacts(end+1,:) = [locsR(i),pksR(i)];
    end
end
for i = 2:length(pksL)
    btw_pk = fsrD(locsL(i-1):locsL(i),Mfsr('Lheel'));
    if any(btw_pk < 1)
        impacts(end+1,:) = [locsL(i),pksL(i)];
    end
end
impacts = sortrows(impacts);

% visually check
figure; hold on
plot(fsrD(:,Mfsr('Lheel')))
plot(fsrD(:,Mfsr('Rheel')))
plot(impacts(:,1),impacts(:,2), 'ko','MarkerSize',10)

%% extract pcb segments and mocap locations from fsr heel impacts

impactT = fsrT(impacts(:,1));
num_sensors = 4;

% pcb extract
window_width = Fs*0.4; % shaking lasts < 0.2s
offset = 0.2; % start window 1/4 of window behind the fsr start
arrival_idx = zeros(length(impactT),num_sensors);
peak_idx = zeros(length(impactT),num_sensors);
peak_mag = zeros(length(impactT),num_sensors);
for i = 1:length(impactT)
    pcbi = findTindex(impactT(i),pcbT);
    starti = pcbi - offset*window_width;
    for j = 1:num_sensors
        endi = min(starti+window_width, length(filt_pcbD(:,1)));
        window = filt_pcbD(starti:endi,j);
        if any(window > 0)
            arrival_idx(i,j) = aic_pick(window, 'to_peak')+starti;
            [mag,idx] = max(abs(window));
            peak_mag(i,j) = mag;
            peak_idx(i,j) = idx + starti;
        else
            arrival_idx(i,j) = NaN;
            peak_idx(i,j) = NaN;
            peak_mag(i,j) = NaN;
        end
    end
end

% mocap extract
footstep_time = 0.3; % length of time foot is reliably still
coordinatesR = ones(length(impactT),2);
coordinatesL = ones(length(impactT),2);
whichfoot = zeros(length(impactT),1);
extracted_pts_R = ones(length(impactT),1);
extracted_pts_L = ones(length(impactT),1);
rightft = 1;
leftft = 2;
avg_seg = 10; % take the average of 10 datapoints
for i = 1:length(impactT)
    starti = findTindex(impactT(i),mocapT);
    endi = findTindex(impactT(i)+footstep_time,mocapT);
    if isempty(endi)
        endi = length(mocapT);
    end
    segL = mocapL(starti:endi,:);
    segR = mocapR(starti:endi,:);
    if std(segL(:,1)) > std(segR(:,1)) % left seg has more variation
        % right leg is still
        prev_avg = mean(segR(1:1+avg_seg,1)); % loop until x-value is steady
        for j = avg_seg:(length(segR)-avg_seg)
            curr_avg = mean(segR(j:j+avg_seg,1));
            if abs(prev_avg-curr_avg) < 50 % movement is < 50 mm
                coordinatesR(i,1) = segR(j,1);
                coordinatesR(i,2) = segR(j,3); % store x,z coordinates
                whichfoot(i,1) = rightft;
                extracted_pts_R(i,1) = starti + j;
                break;
            end
            j = j+avg_seg;
        end
    else
        % left leg is still
        prev_avg = mean(segL(1:1+avg_seg,1)); % loop until x-value is steady
        for j = avg_seg:(length(segL)-avg_seg)
            curr_avg = mean(segL(j:j+avg_seg,1));
            if abs(prev_avg-curr_avg) < 50 % movement is < 50 mm
                coordinatesL(i,1) = segL(j,1);
                coordinatesL(i,2) = segL(j,3); % store x,z coordinates
                whichfoot(i,1) = leftft;
                extracted_pts_L(i,1) = starti + j;
                break;
            end
            j = j+avg_seg;
        end
    end
end

% visualize
peak_idx(isnan(peak_idx)) = 1;
peak_mag(isnan(peak_mag)) = 1;
arrival_idx(isnan(arrival_idx)) = 1;
figure;
for i=1:num_sensors
    subplot(4,1,i)
    plot(pcbT,filt_pcbD(:,i))
    hold on
    plot(pcbT(peak_idx(:,i)),peak_mag(:,i),'r.','MarkerSize',10)
    plot(pcbT(arrival_idx(:,i)),0,'k.','MarkerSize',10)
end

figure;
subplot(2,1,1)
plot(mocapT,mocapR(:,1))
hold on
plot(mocapT(extracted_pts_R), coordinatesR(:,1),'r.','MarkerSize',10)
subplot(2,1,2)
plot(mocapT,mocapR(:,3))
hold on
plot(mocapT(extracted_pts_R), coordinatesR(:,2),'r.','MarkerSize',10)
title('Right foot mocap')

figure;
subplot(2,1,1)
plot(mocapT,mocapL(:,1))
hold on
plot(mocapT(extracted_pts_L), coordinatesL(:,1),'r.','MarkerSize',10)
subplot(2,1,2)
plot(mocapT,mocapL(:,3))
hold on
plot(mocapT(extracted_pts_L), coordinatesL(:,2),'r.','MarkerSize',10)
title('Left foot mocap')

%% probably trash everything under this: was used to make a small plot for mini group
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




