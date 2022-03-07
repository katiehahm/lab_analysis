data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\ProcessedData\both_straight2';
load(string(data_root_katie))

figure;
plot(pcbTime, filt_pcbD(:,1))
hold on
plot(fsrTime(impacts(:,1)), 0, 'rx')

% clip 33.24 s to 38.8302
pcb_clip = filt_pcbD(425474:497025,1);
time_clip = pcbTime(425474:497025);
Fs = 12800;
% helperPlotSpectrogram(pcb_clip,time_clip,Fs,200)
% 
% spectrogram(pcb_clip,256,250,[],Fs,'yaxis')
% cwt(pcb_clip,Fs)
% 
% helperPlotScalogram(pcb_clip,Fs)
% stft(pcb_clip,Fs)

%%
figure;
spectrogram(pcb_clip,256,250,[],Fs,'yaxis');
%%
figure;
plot(pcbTime, rawPcb(:,3))
hold on
plot(fsrTime(impacts(:,1)), 0, 'rx')

%% 3/6/22 assuming that we can detect all arrival times, try to parse & label

impact_starts = fsrTime(impacts(:,1));
impact_starts = sort(impact_starts); % time of start of impacts

impact_ends = get_impact_window(impact_starts, filt_pcbD, pcbTime); % returns sensorN x impactN matrix of impact end times

figure; plot(pcbTime, filt_pcbD(:,1))
hold on
plot(impact_starts,0,'rx','MarkerSize',10)
hold on
plot(impact_ends(:,1),0,'cx','MarkerSize',10)

%% 3/7/22 xcorr to label impacts A or B

% looks at the next 3 impacts to label them
real_ID = zeros(length(impact_starts),1);
person1_idx = find(impacts(:,4) == 1 || impacts(:,4) == 2);
real_ID(person1_idx) = 1; % label person 1 as 1
person2_idx = find(impacts(:,4) == 3 || impacts(:,4) == 4);
real_ID(person2_idx) = 2; % label person 2 as 2

xcorr_estimates = zeros(length(impact_starts),






















