% 3/10/22 uses crc data to see correlation between grf and ta

clear all
close all
fileroot = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\ExperimentCRC\';
take = '5'; 
filepath = [fileroot, 'Take',take];
% load datasets
load([filepath, '_delsys'])
load([filepath,])

tap_fsrID = 1;

% Fs
Fs_fsr1 = Fs(1);
Fs_fsr = Fs(2);
Fs_acc = Fs(3);

% process struct from force plates
fp_struct = Katie_03180005; % change this value #######
sCell = struct2cell(fp_struct.Force);
Fs_fp = cell2mat(sCell(5,1,1)); % Fs
fpN = cell2mat(sCell(4,1,1)); % number of frames
forceX = zeros(4,fpN);
forceY = zeros(4,fpN);
forceZ = zeros(4,fpN);
for i = 1:4
    forces = cell2mat(sCell(6,1,i));
    forceX(i,:) = forces(1,:);
    forceY(i,:) = forces(2,:);
    forceZ(i,:) = forces(3,:);
end

figure;
plot(forceZ(2,:)) % taps are on plate 2
figure;
plot(Data(tap_fsrID,:))

%% synchronize

fsr_peaks = [16254,16453,16672,16826];
fp_peaks = [19462,19959,20424,20846];

for i = 2:length(fp_peaks)
    (fsr_peaks(i) - fsr_peaks(i-1))/Fs_fsr1
    (fp_peaks(i) - fp_peaks(i-1))/Fs_fp
end

%% cut off at peak
peak_num = 3;

forceX = forceX(:,fp_peaks(peak_num)+Fs_fp:end);
forceY = forceY(:,fp_peaks(peak_num)+Fs_fp:end);
forceZ = forceZ(:,fp_peaks(peak_num)+Fs_fp:end);
fpTime = linspace(0,fpN/Fs_fp,fpN);
fpTime = fpTime(fp_peaks(peak_num)+Fs_fp:end) - fpTime(fp_peaks(peak_num)+Fs_fp);

trigger_num = fsr_peaks(peak_num);
fsr_data_num = trigger_num*Fs_fsr/Fs_fsr1;
accel_data_num = trigger_num*Fs_acc/Fs_fsr1;

fsrData = Data(2,round(fsr_data_num+Fs_fsr):end); % left foot
fsrData = [fsrData; Data(6,round(fsr_data_num+Fs_fsr):end)]; % right foot
fsrTime = Time(2,:);
fsrTime = fsrTime(round(fsr_data_num + Fs_fsr):end)-fsrTime(round(fsr_data_num + Fs_fsr));

accData = Data(3:5,round(accel_data_num+Fs_acc):end); % left foot
accData = [accData; Data(7:9,round(accel_data_num+Fs_acc):end)]; % left foot
accTime = Time(3,:);
accTime = accTime(round(accel_data_num + Fs_acc):end)-accTime(round(accel_data_num + Fs_acc));

triggerData = Data(1,round(trigger_num + Fs_fsr1):end);

% just make sure it's cut properly
figure; plot(triggerData)
figure; plot(forceZ(2,:))

% for next section
figure; plot(sum(forceZ,1))

%% find start and end time

start_idx = 12401; % index from fp
last_idx = 333289; %index from fp
start_time = start_idx/Fs_fp;
last_time = last_idx/Fs_fp;

forceX = forceX(:,start_idx:last_idx);
forceY = forceY(:,start_idx:last_idx);
forceZ = forceZ(:,start_idx:last_idx);
fpTime = fpTime(start_idx:last_idx);

start_fsr_idx = findTindex(start_time, fsrTime);
last_fsr_idx = findTindex(last_time, fsrTime);
start_acc_idx = findTindex(start_time, accTime);
last_acc_idx = findTindex(last_time, accTime);
fsrData = fsrData(:,start_fsr_idx:last_fsr_idx);
fsrTime = fsrTime(start_fsr_idx:last_fsr_idx);
accData = accData(:,start_acc_idx:last_acc_idx);
accTime = accTime(start_acc_idx:last_acc_idx);

figure;
subplot(3,1,1)
plot(fpTime, forceZ(1,:))
hold on
plot(fpTime, forceZ(2,:))
plot(fpTime, forceZ(3,:))
plot(fpTime, forceZ(4,:))
subplot(3,1,2)
plot(fsrTime, fsrData(1,:),'r')
hold on
plot(fsrTime, fsrData(2,:),'b')
subplot(3,1,3)
plot(accTime, accData(2,:),'r')
hold on
plot(accTime, accData(5,:),'b')

%% save data
processedfilepath = [fileroot,'ProcessedData\Take',take];

save(processedfilepath,'accData','accTime','Fs_acc','Data',...
    'forceX','forceY','forceZ','fpTime','Fs_fp','fpN',...
    'fsrData','fsrTime','Fs_fsr','Time','Katie_03180005')
disp(append("Saved as ", processedfilepath))

%% find impact times & magnitudes

% use GRF impacts to find peaks, then use those timestamps to find fsr time
% to get accel peaks

close all
fp_peaks = [0,0,0];
for i = 1:4 % 4 force plates
    figure; findpeaks(forceZ(i,:),'MinPeakProminence',300,'MinPeakDistance',3*Fs_fp)
    [pks,locs] = findpeaks(forceZ(i,:),'MinPeakProminence',300,'MinPeakDistance',3*Fs_fp);
    identifier = zeros(length(locs),1);
    identifier(:) = i;
    add_matrix = [locs.', pks.',identifier];
    fp_peaks = [fp_peaks; add_matrix]; % each row is new impact, each row is (loc, pk)
end

fp_peaks(1,:) = []; % delete first row from initialization
fp_peaks = sortrows(fp_peaks,1);
impactN = length(fp_peaks);

% save(processedfilepath,'fp_peaks','impactN','-append')



%% use these peaks to find fsr timing and accel peaks

clear all
load('C:\Users\katie\Dropbox (MIT)\Lab\Analysis\ExperimentCRC\ProcessedData\Take3.mat')

window = 0.25; % 0.25 around the fp impact to find impact start in fsr

impacts = zeros(impactN,9); % refer to notebook

for i = 1:impactN
    curr_time = fp_peaks(i,1)/Fs_fp;
    start_idx = round((curr_time - window*0.5)*Fs_fsr);
    acc_start_idx = round(start_idx*Fs_acc/Fs_fsr);
    left_window = fsrData(1,start_idx:start_idx + round(window*Fs_fsr));
    right_window = fsrData(2,start_idx:start_idx + round(window*Fs_fsr));
    if max(left_window) > max(right_window) % impact happened on the left foot
        [peak, idx] = max(left_window);
        labeler = 1;
    else % right foot
        [peak, idx] = max(right_window);
        labeler = 2;
    end
    acc_windowX = accData((labeler-1)*3 + 1,acc_start_idx:acc_start_idx + round(window*Fs_acc)); % tibial x-axis
    acc_windowY = accData((labeler-1)*3 + 2,acc_start_idx:acc_start_idx + round(window*Fs_acc)); % tibial y-axis
    acc_windowZ = accData((labeler-1)*3 + 3,acc_start_idx:acc_start_idx + round(window*Fs_acc)); % tibial z-axis
    [accpeakX, accidxX] = max(abs(acc_windowX));
    [accpeakY, accidxY] = max(abs(acc_windowY));
    [accpeakZ, accidxZ] = max(abs(acc_windowZ));
    impacts(i,:) = [idx + start_idx, peak, accidxX + start_idx, accpeakX, accidxY + start_idx, accpeakY, accidxZ + start_idx, accpeakZ, labeler];
end

under8idx = find(impacts(:,6) < 7.9); % this only considers the y axis
X = impacts(under8idx,4)+impacts(under8idx,6)+impacts(under8idx,8);
% X = impacts(under8idx,6); % just y accel
y = fp_peaks(under8idx,2);
figure;
plot(y,X,'b.')
xlabel('Force plate')
ylabel('Tibial Acceleration')

mdl = fitlm(X,y);
anova(mdl,'summary')
figure; plot(mdl)

mdl
mdl.RMSE/mean(y)


%% finding the first peak of fp for each impact 3/19/22

% use fp_peaks to find the first peak (FP data generally has 2 peaks)
% get window (by using threshold < 50)
% cut it in half (they all seem to have two humps)
% get the max value

fp_window_size = 1300; % data samples
threshold = 50;
fp_first_peaks = [0,0,0];
for i = 1:impactN
    fp_label = fp_peaks(i,3);
    fp_idx = fp_peaks(i,1);
    start_idx = max(fp_idx - fp_window_size,1);
    fp_window = forceZ(fp_label,start_idx:min(fp_idx+fp_window_size,length(forceZ)));
    include_idx = find(fp_window > threshold);
    fp_window = fp_window(include_idx);
    fp_window = fp_window(1:round(length(fp_window)/2));
    [peak, idx] = max(fp_window);
    add_matrix = [start_idx + include_idx(1) - 1 + idx - 1, peak, fp_label];
    fp_first_peaks = [fp_first_peaks; add_matrix];
end

fp_first_peaks(1,:) = []; % delete first row from initialization

for i = 1:4
    figure;
    plot(forceZ(i,:))
    hold on
    first_idx = find(fp_first_peaks(:,3) == i);
    plot(fp_first_peaks(first_idx,1),fp_first_peaks(first_idx,2),'rx')
end

% save(processedfilepath,'fp_peaks','impactN','-append')

%% use fp_first_peaks to find fsr timing and accel peaks

window = 0.25; % 0.25 around the fp impact to find impact start in fsr

impacts = zeros(impactN,9); % refer to notebook

for i = 1:impactN
    curr_time = fp_first_peaks(i,1)/Fs_fp;
    start_idx = round((curr_time - window*0.5)*Fs_fsr);
    acc_start_idx = round(start_idx*Fs_acc/Fs_fsr);
    left_window = fsrData(1,start_idx:start_idx + round(window*Fs_fsr));
    right_window = fsrData(2,start_idx:start_idx + round(window*Fs_fsr));
    if max(left_window) > max(right_window) % impact happened on the left foot
        [peak, idx] = max(left_window);
        labeler = 1;
    else % right foot
        [peak, idx] = max(right_window);
        labeler = 2;
    end
    acc_windowX = accData((labeler-1)*3 + 1,acc_start_idx:acc_start_idx + round(window*Fs_acc)); % tibial x-axis
    acc_windowY = accData((labeler-1)*3 + 2,acc_start_idx:acc_start_idx + round(window*Fs_acc)); % tibial y-axis
    acc_windowZ = accData((labeler-1)*3 + 3,acc_start_idx:acc_start_idx + round(window*Fs_acc)); % tibial z-axis
    [accpeakX, accidxX] = max(abs(acc_windowX));
    [accpeakY, accidxY] = max(abs(acc_windowY));
    [accpeakZ, accidxZ] = max(abs(acc_windowZ));
    impacts(i,:) = [idx + start_idx, peak, accidxX + start_idx, accpeakX, accidxY + start_idx, accpeakY, accidxZ + start_idx, accpeakZ, labeler];
end

under8idx = find(impacts(:,6) < 7.9); % this only considers the y axis
X = impacts(under8idx,4)+impacts(under8idx,6)+impacts(under8idx,8);
% X = impacts(under8idx,6); % just y accel
y = fp_first_peaks(under8idx,2);
figure;
plot(y,X,'b.')
xlabel('Force plate')
ylabel('Tibial Acceleration')

mdl = fitlm(X,y);
anova(mdl,'summary')
figure; plot(mdl)

mdl
mdl.RMSE/mean(y)




