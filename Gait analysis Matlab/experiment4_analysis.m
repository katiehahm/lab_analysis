% this file was used to experiment with footfall detection
% does not contain any reusable code

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Praneeth experiment 3_11_22\ProcessedData\both_towards1';
load(string(data_root_katie))

figure;
plot(pcbTime, filt_pcbD(:,1))
hold on
plot(fsrTime(impacts(:,1)), 0, 'rx')

%% make clip
start_time = 41.8;
end_time = 55.43;
clip_start = findTindex(start_time,pcbTime);
clip_end = findTindex(end_time,pcbTime);
pcb_clip = filt_pcbD(clip_start:clip_end,1);
time_clip = pcbTime(clip_start:clip_end);
impact_idx = find(fsrTime(impacts(:,1)) < end_time & fsrTime(impacts(:,1)) > start_time);
clip_impact_times = fsrTime(impacts(impact_idx,1)) - start_time;

figure;
plot(time_clip - time_clip(1), pcb_clip)
hold on
plot(clip_impact_times,0,'rx')
title('Signal clip example')
% helperPlotSpectrogram(pcb_clip,time_clip,Fs,200)
% 
% spectrogram(pcb_clip,256,250,[],Fs,'yaxis')
% cwt(pcb_clip,Fs)
% 
% helperPlotScalogram(pcb_clip,Fs)
% stft(pcb_clip,Fs)

%%
figure;
spectrogram(filt_pcbD(clip_start:clip_end,5),256,250,[],Fs_pcb,'yaxis');
hold on
plot(clip_impact_times,0,'rx')
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

% [impactN, sensorN] = size(arrival_idx);
% % looks at the next 3 impacts to label them
% real_ID = zeros(impactN,1);
% person1_idx = find(impacts(:,4) == 1 || impacts(:,4) == 2);
% real_ID(person1_idx) = 1; % label person 1 as 1
% person2_idx = find(impacts(:,4) == 3 || impacts(:,4) == 4);
% real_ID(person2_idx) = 2; % label person 2 as 2
% 
% xcorr_estimates = zeros(impactN,sensorN);
% xcorr_estimates(1,:) = 1; % always label 1st impact as person 1
% 
% for s = 1:sensorN
%     for i = 1:impactN
%         if (impactN-i) == 3 % 4th to last
%             
%         elseif (impactN-i) == 2
%             
%         else
%             curr_ID = xcorr_estimates(i,s);
%             curr_window = filt_pcbD(findTindex(impact_starts(i),pcbTime):findTindex(impact_ends(i,s),pcbTime));
%             first
%         end
%     end
% end

%% 3/7/22 xcorr doesn't work try fft
% 
% % example data, clip 1 & 3, clip 2 & 4 are same person (v obvious)
% clip1 = filt_pcbD(findTindex(11.3201,pcbTime):findTindex(11.4409,pcbTime));
% clip2 = filt_pcbD(findTindex(11.5564,pcbTime):findTindex(11.6913,pcbTime));
% clip3 = filt_pcbD(findTindex(11.7589,pcbTime):findTindex(11.8772,pcbTime));
% clip4 = filt_pcbD(findTindex(12.0998,pcbTime):findTindex(12.2246,pcbTime));
% 
% % HERE look at fft example to plot them
% fft_clip1 = fft(clip1);
% L_clip1 = length(clip1);
% P2_clip1 = abs(fft_clip1/L_clip1);
% P1_clip1 = P2_clip1(1:L_clip1/2+1);
% P1_clip1(2:end-1) = 2*P1_clip1(2:end-1);
% 
% f_clip1 = Fs*(0:(L_clip1/2))/L_clip1;
% figure; plot(f_clip1,P1_clip1)
% title('Single-Sided Amplitude Spectrum of clip 1')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% 
% fft_clip2 = fft(clip2);
% L_clip2 = length(clip2);
% P2_clip2 = abs(fft_clip2/L_clip2);
% P1_clip2 = P2_clip2(1:L_clip2/2+1);
% P1_clip2(2:end-1) = 2*P1_clip2(2:end-1);
% 
% f_clip2 = Fs*(0:(L_clip2/2))/L_clip2;
% figure; plot(f_clip2,P1_clip2)
% title('Single-Sided Amplitude Spectrum of clip 2')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% 
% fft_clip3 = fft(clip3);
% L_clip3 = length(clip3);
% P2_clip3 = abs(fft_clip3/L_clip3);
% P1_clip3 = P2_clip3(1:L_clip3/2+1);
% P1_clip3(2:end-1) = 2*P1_clip3(2:end-1);
% 
% f_clip3 = Fs*(0:(L_clip3/2))/L_clip3;
% figure; plot(f_clip3,P1_clip3)
% title('Single-Sided Amplitude Spectrum of clip 3')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% 
% fft_clip4 = fft(clip4);
% L_clip4 = length(clip4);
% P2_clip4 = abs(fft_clip4/L_clip4);
% P1_clip4 = P2_clip4(1:L_clip4/2+1);
% P1_clip4(2:end-1) = 2*P1_clip4(2:end-1);
% 
% f_clip4 = Fs*(0:(L_clip4/2))/L_clip4;
% figure; plot(f_clip4,P1_clip4)
% title('Single-Sided Amplitude Spectrum of clip 4')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% 
% % find fft up to 2000 Hz (signal 0 after that), interpolate, then multiply
% % fft to compare similarity
% commonF = linspace(0,2000,1000); % increments of 5 hz
% cutoff1 = findTindex(2000,f_clip1);
% P1_clip1_common = interp1(f_clip1(1:cutoff1),P1_clip1(1:cutoff1),commonF);
% figure; plot(commonF, P1_clip1_common)
% title('Single-Sided Amplitude Spectrum of clip 1 Interpolated')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% 
% cutoff2 = findTindex(2000,f_clip2);
% P1_clip2_common = interp1(f_clip2(1:cutoff2),P1_clip2(1:cutoff2),commonF);
% figure; plot(commonF, P1_clip2_common)
% title('Single-Sided Amplitude Spectrum of clip 2 Interpolated')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% 
% cutoff3 = findTindex(2000,f_clip3);
% P1_clip3_common = interp1(f_clip3(1:cutoff3),P1_clip3(1:cutoff3),commonF);
% figure; plot(commonF, P1_clip3_common)
% title('Single-Sided Amplitude Spectrum of clip 3 Interpolated')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% 
% cutoff4 = findTindex(2000,f_clip4);
% P1_clip4_common = interp1(f_clip4(1:cutoff4),P1_clip4(1:cutoff4),commonF);
% figure; plot(commonF, P1_clip4_common)
% title('Single-Sided Amplitude Spectrum of clip 4 Interpolated')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% 
% %% fft cont.
% 
% P1_norm1 = P1_clip1_common./max(P1_clip1_common);
% P1_norm2 = P1_clip2_common./max(P1_clip2_common);
% P1_norm3 = P1_clip3_common./max(P1_clip3_common);
% P1_norm4 = P1_clip4_common./max(P1_clip4_common);
% 
% clip1and2 = sum(P1_norm1.*P1_norm2)
% clip1and3 = sum(P1_norm1.*P1_norm3)
% clip1and4 = sum(P1_norm1.*P1_norm4)
% clip2and3 = sum(P1_norm2.*P1_norm3)
% clip2and4 = sum(P1_norm2.*P1_norm4)
% clip3and4 = sum(P1_norm3.*P1_norm4)
% 
% r12 = xcorr(P1_norm1,P1_norm2);
% r13 = xcorr(P1_norm1,P1_norm3);
% r14 = xcorr(P1_norm1,P1_norm4);
% r23 = xcorr(P1_norm2,P1_norm3);
% r24 = xcorr(P1_norm2,P1_norm4);
% r34 = xcorr(P1_norm3,P1_norm4);
% max(r12)
% max(r13)
% max(r14)
% max(r23)
% max(r24)
% max(r34)
% 
% clip1and2 = sum(P1_clip1_common.*P1_clip2_common)
% clip1and3 = sum(P1_clip1_common.*P1_clip3_common)
% clip1and4 = sum(P1_clip1_common.*P1_clip4_common)
% clip2and3 = sum(P1_clip2_common.*P1_clip3_common)
% clip2and4 = sum(P1_clip2_common.*P1_clip4_common)
% clip3and4 = sum(P1_clip3_common.*P1_clip4_common)
% 
% energy1 = sum(P1_norm1)
% energy2 = sum(P1_norm2)
% energy3 = sum(P1_norm3)
% energy4 = sum(P1_norm4)

% %% new impact ends (shorter window if next impact is within window length & peak is smaller than next impact peak 3/10/22
% 
% new_impact_ends = zeros(impactN,sensorN);
% for s = 1:sensorN
%     for i = 1:impactN-1
%         if impact_starts(i+1) < impact_ends(i,s)
%             
%         end
%     end
% end


%% fft cont. using energy to differentiate 3/8/22

Fs = 12800;
[impactN, sensorN] = size(arrival_idx);
% looks at the next 3 impacts to label them
real_ID = zeros(impactN,1);
person1_idx = find(impacts(:,4) == 1 | impacts(:,4) == 2);
real_ID(person1_idx) = 1; % label person 1 as 1
person2_idx = find(impacts(:,4) == 3 | impacts(:,4) == 4);
real_ID(person2_idx) = 2; % label person 2 as 2

ID_estimates = zeros(impactN,sensorN);
ID_estimates(1,:) = 2; % always label 1st impact as person 1
ID_estimates(2,:) = 1; % always label 2nd impact as person 2
ID_confidence = zeros(impactN,sensorN);

for s = 1:sensorN
    for i = 1:impactN
        curr_ID = ID_estimates(i,s);
        curr_clip = filt_pcbD(findTindex(impact_starts(i),pcbTime):findTindex(impact_ends(i,s),pcbTime),s);
        first_clip = [];
        sec_clip = [];
        third_clip = [];
        commonF = linspace(0,2000,1000);
        
        if i+1 <= impactN && ID_estimates(i+1,s) == 0
            first_clip = filt_pcbD(findTindex(impact_starts(i+1),pcbTime):findTindex(impact_ends(i+1,s),pcbTime),s);
        end
        if i+2 <= impactN && ID_estimates(i+2,s) == 0
            sec_clip = filt_pcbD(findTindex(impact_starts(i+2),pcbTime):findTindex(impact_ends(i+2,s),pcbTime),s);
        end
        if i+3 <= impactN && ID_estimates(i+3,s) == 0
            third_clip = filt_pcbD(findTindex(impact_starts(i+3),pcbTime):findTindex(impact_ends(i+3,s),pcbTime),s);
        end
        if ~isempty([first_clip;sec_clip;third_clip])
            clip_energy = fft_energy(curr_clip,first_clip,sec_clip,third_clip,commonF,Fs);
            energy_ratios = [];
            for k = 2:length(clip_energy)
                if clip_energy(k) ~= 0
                    energy_ratios(end+1) = clip_energy(k)/clip_energy(1);
                else
                    energy_ratios(end+1) = NaN;
                end
            end
            [val,idx] = min(abs(energy_ratios - 1)); % get the idx with ratio closest to 1
            ID_estimates(i+idx,s) = curr_ID;
            % confidence calculates how different it is from max value
            ID_confidence(i+idx,s) = abs(val - max(abs(energy_ratios - 1)));               
        end
    end
end

final_estimates = zeros(impactN,1);
% make sensors vote for the most likely ID
% take vote of the lowest uncertainty
for i = 1:impactN
    curr_estimates = ID_estimates(i,:);
    curr_confidence = ID_confidence(i,:);
    
%     final_estimates(i) = mode(curr_estimates); % vote most often 

%     [~,max_id] = max(curr_confidence); % using max confidence
%     final_estimates(i) = curr_estimates(max_id);
    
    vote1 = find(curr_estimates == 1); % sum of confidences
    vote2 = find(curr_estimates == 2);
    vote1confidence = sum(curr_confidence(vote1));
    vote2confidence = sum(curr_confidence(vote2));
    if vote1confidence > vote2confidence | isempty(vote2)
        
        final_estimates(i) = 1;
    else
        final_estimates(i) = 2;
    end
end

figure;
plot(person1_idx,zeros(length(person1_idx),1),'rx')
hold on
plot(person2_idx,zeros(length(person2_idx),1),'bo')
est_person1_idx = find(final_estimates == 1);
est_person2_idx = find(final_estimates == 2);
plot(est_person1_idx,ones(length(est_person1_idx),1),'rx')
plot(est_person2_idx,ones(length(est_person2_idx),1),'bo')
legend('Real person 1','Real person 2','Estimated person 1','Estimated person 2')
ylim([-1 2])

figure;
results = real_ID - final_estimates;
plot(results,'.')
title('Real ID - estimated ID')
ylim([-2 2])
correct_idx = find(results == 0);
accuracy = length(correct_idx)/length(results)

%% similar approach to above (fft energy) but with just magnitude 3/15/22

Fs = 12800;
[impactN, sensorN] = size(arrival_idx);
% looks at the next 3 impacts to label them
real_ID = zeros(impactN,1);
person1_idx = find(impacts(:,4) == 1 | impacts(:,4) == 2);
real_ID(person1_idx) = 1; % label person 1 as 1
person2_idx = find(impacts(:,4) == 3 | impacts(:,4) == 4);
real_ID(person2_idx) = 2; % label person 2 as 2

ID_estimates = zeros(impactN,sensorN);
ID_estimates(1,:) = 2; % always label 1st impact as person 1
ID_estimates(2,:) = 1; % always label 2nd impact as person 2
% ID_estimates(3,:) = 2;
ID_confidence = zeros(impactN,sensorN);

% for i = 1:length(impact_starts)-1
%     if impact_starts(i) == impact_starts(i+1)
%         impact_starts(i) = impact_starts(i) - 5/Fs_pcb;
%     end
% end
% 
% impact_ends = old_impact_ends;
% for s = 1:sensorN
%     for i = 1:length(impact_starts)-1
%         if impact_ends(i,s) > impact_starts(i+1)
%             first_max = max(filt_pcbD(findTindex(impact_starts(i),pcbTime):findTindex(impact_starts(i+1),pcbTime),s));
%             next_max = max(filt_pcbD(findTindex(impact_starts(i+1),pcbTime):findTindex(impact_ends(i+1,s),pcbTime),s));
%             if first_max < next_max % nearly overlapping two impacts have the first one smaller
%                 impact_ends(i,s) = impact_starts(i+1);
%             end
%         end
%     end
% end



for s = 1:sensorN
    for i = 1:impactN
        curr_ID = ID_estimates(i,s);
        curr_clip = filt_pcbD(findTindex(impact_starts(i),pcbTime):findTindex(impact_ends(i,s),pcbTime),s);
        first_clip = [];
        sec_clip = [];
        third_clip = [];
        commonF = linspace(0,2000,1000);
        
        if i+1 <= impactN && ID_estimates(i+1,s) == 0
            first_clip = filt_pcbD(findTindex(impact_starts(i+1),pcbTime):findTindex(impact_ends(i+1,s),pcbTime),s);
        end
        if i+2 <= impactN && ID_estimates(i+2,s) == 0
            sec_clip = filt_pcbD(findTindex(impact_starts(i+2),pcbTime):findTindex(impact_ends(i+2,s),pcbTime),s);
        end
        if i+3 <= impactN && ID_estimates(i+3,s) == 0
            third_clip = filt_pcbD(findTindex(impact_starts(i+3),pcbTime):findTindex(impact_ends(i+3,s),pcbTime),s);
        end
        
        if ~isempty([first_clip;sec_clip;third_clip])
            max_curr = max(abs(curr_clip));
            if ~isempty(first_clip)
                max_first = max(abs(first_clip));
            else
                max_first = NaN;
            end
            if ~isempty(sec_clip)
                max_sec = max(abs(sec_clip));
            else
                max_sec = NaN;
            end
            if ~isempty(third_clip)
                max_third = max(abs(third_clip));
            else
                max_third = NaN;
            end
            max_amplitude = [max_curr, max_first, max_sec, max_third];
            amp_ratios = [];
            for k = 2:length(max_amplitude)
                if ~isempty(max_amplitude(k))
                    amp_ratios(end+1) = max_amplitude(k)/max_amplitude(1);
                else
                    amp_ratios(end+1) = NaN;
                end
            end
            [val,idx] = min(abs(amp_ratios - 1)); % get the idx with ratio closest to 1
            ID_estimates(i+idx,s) = curr_ID;
            % confidence calculates how different it is from max value
            ID_confidence(i+idx,s) = abs(val - max(abs(amp_ratios - 1)));
        end
    end
end

final_estimates = zeros(impactN,1);
% make sensors vote for the most likely ID
% take vote of the lowest uncertainty
for i = 1:impactN
    curr_estimates = ID_estimates(i,:);
    curr_confidence = ID_confidence(i,:);
    
%     final_estimates(i) = mode(curr_estimates); % vote most often 

%     [~,max_id] = max(curr_confidence); % using confidence
%     final_estimates(i) = curr_estimates(max_id);
    
    vote1 = find(curr_estimates == 1); % sum of confidences
    vote2 = find(curr_estimates == 2);
    vote1confidence = sum(curr_confidence(vote1));
    vote2confidence = sum(curr_confidence(vote2));
    if vote1confidence > vote2confidence
        final_estimates(i) = 1;
    else
        final_estimates(i) = 2;
    end
end

figure;
plot(person1_idx,zeros(length(person1_idx),1),'rx')
hold on
plot(person2_idx,zeros(length(person2_idx),1),'bo')
est_person1_idx = find(final_estimates == 1);
est_person2_idx = find(final_estimates == 2);
plot(est_person1_idx,ones(length(est_person1_idx),1),'rx')
plot(est_person2_idx,ones(length(est_person2_idx),1),'bo')
legend('Real person 1','Real person 2','Estimated person 1','Estimated person 2')
ylim([-1 2])

figure;
results = real_ID - final_estimates;
plot(results,'.')
title('Real ID - estimated ID')
ylim([-2 2])
correct_idx = find(results == 0);
accuracy = length(correct_idx)/length(results)

%% severe bandpass filter 3/10/22 (didnt see a big difference)

bp_pcbD = bandpass(filt_pcbD,[150 200],Fs);
figure;
subplot(2,1,1)
plot(pcbTime, filt_pcbD(:,1))
hold on
plot(impact_starts,0,'rx','MarkerSize',12)
subplot(2,1,2)
plot(pcbTime, bp_pcbD(:,1))
hold on
plot(impact_starts,0,'rx','MarkerSize',12)



%% fft cont. using energy with set window length to differentiate 3/17/22

Fs = 12800;
[impactN, sensorN] = size(arrival_idx);
% looks at the next 3 impacts to label them
real_ID = zeros(impactN,1);
person1_idx = find(impacts(:,4) == 1 | impacts(:,4) == 2);
real_ID(person1_idx) = 1; % label person 1 as 1
person2_idx = find(impacts(:,4) == 3 | impacts(:,4) == 4);
real_ID(person2_idx) = 2; % label person 2 as 2

ID_estimates = zeros(impactN,sensorN);
ID_estimates(1,:) = 2; % always label 1st impact as person 1
ID_estimates(2,:) = 1; % always label 2nd impact as person 2
ID_confidence = zeros(impactN,sensorN);

set_impact_ends = impact_starts + 0.25;

for s = 1:sensorN
    for i = 1:impactN
        curr_ID = ID_estimates(i,s);
        curr_clip = filt_pcbD(findTindex(impact_starts(i),pcbTime):findTindex(set_impact_ends(i),pcbTime),s);
        first_clip = [];
        sec_clip = [];
        third_clip = [];
        commonF = linspace(0,2000,1000);
        
        if i+1 <= impactN && ID_estimates(i+1,s) == 0
            first_clip = filt_pcbD(findTindex(impact_starts(i+1),pcbTime):findTindex(set_impact_ends(i+1),pcbTime),s);
        end
        if i+2 <= impactN && ID_estimates(i+2,s) == 0
            sec_clip = filt_pcbD(findTindex(impact_starts(i+2),pcbTime):findTindex(set_impact_ends(i+2),pcbTime),s);
        end
        if i+3 <= impactN && ID_estimates(i+3,s) == 0
            third_clip = filt_pcbD(findTindex(impact_starts(i+3),pcbTime):findTindex(set_impact_ends(i+3),pcbTime),s);
        end
        if ~isempty([first_clip;sec_clip;third_clip])
            clip_energy = fft_energy(curr_clip,first_clip,sec_clip,third_clip,commonF,Fs);
            energy_ratios = [];
            for k = 2:length(clip_energy)
                if clip_energy(k) ~= 0
                    energy_ratios(end+1) = clip_energy(k)/clip_energy(1);
                else
                    energy_ratios(end+1) = NaN;
                end
            end
            [val,idx] = min(abs(energy_ratios - 1)); % get the idx with ratio closest to 1
            ID_estimates(i+idx,s) = curr_ID;
            % confidence calculates how different it is from max value
            ID_confidence(i+idx,s) = abs(val - max(abs(energy_ratios - 1)));               
        end
    end
end

final_estimates = zeros(impactN,1);
% make sensors vote for the most likely ID
% take vote of the lowest uncertainty
for i = 1:impactN
    curr_estimates = ID_estimates(i,:);
    curr_confidence = ID_confidence(i,:);
    
%     final_estimates(i) = mode(curr_estimates); % vote most often 

%     [~,max_id] = max(curr_confidence); % using max confidence
%     final_estimates(i) = curr_estimates(max_id);
    
    vote1 = find(curr_estimates == 1); % sum of confidences
    vote2 = find(curr_estimates == 2);
    vote1confidence = sum(curr_confidence(vote1));
    vote2confidence = sum(curr_confidence(vote2));
    if vote1confidence > vote2confidence | isempty(vote2)
        
        final_estimates(i) = 1;
    else
        final_estimates(i) = 2;
    end
end

figure;
plot(person1_idx,zeros(length(person1_idx),1),'rx')
hold on
plot(person2_idx,zeros(length(person2_idx),1),'bo')
est_person1_idx = find(final_estimates == 1);
est_person2_idx = find(final_estimates == 2);
plot(est_person1_idx,ones(length(est_person1_idx),1),'rx')
plot(est_person2_idx,ones(length(est_person2_idx),1),'bo')
legend('Real person 1','Real person 2','Estimated person 1','Estimated person 2')
ylim([-1 2])

figure;
results = real_ID - final_estimates;
plot(results,'.')
title('Real ID - estimated ID')
ylim([-2 2])
correct_idx = find(results == 0);
accuracy = length(correct_idx)/length(results)

%% plot the distribution of all step times
% clustering around o-x, x-o, x-x, o-x?

step_times = zeros(length(impact_starts)-1,1);
for i = 1:length(impact_starts)-1
    step = impact_starts(i+1)-impact_starts(i);
    if step < 1
        step_times(i) = step;
    end
end
figure;
histogram(step_times)

%% using step time to ID footsteps 3/20/22

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Praneeth experiment 3_11_22\ProcessedData\both_towards1';
load(string(data_root_katie))

% get the average step time for each person
% and get walk edges
impactN = length(impacts);
walk_edges = zeros(1,impactN);
walk_edges(1) = -1;
walk_edges(end) = 1;
ID_labels = zeros(1,impactN);
ID_o = 0;
ID_o_count = 0;
last_o_idx = 0;
ID_x = 0;
ID_x_count = 0;
last_x_idx = 0;

% this assumes the first two footsteps of every episode is different person
for i = 1:2
    personID = impacts(i,4);
    if personID == 1 | personID == 2
        ID_labels(i) = 1;
        last_o_idx = i;
    else
        ID_labels(i) = 2;
        last_x_idx = i;
    end
end
    

for i = 3:impactN
    if (impacts(i,1) - impacts(i-1,1))/Fs_fsr > 1 % start of new walking segment
        walk_edges(i-1) = 1;
        walk_edges(i) = -1;
        % first step of segment
        personID = impacts(i,4);
        if personID == 1 | personID == 2
            ID_labels(i) = 1;
            last_o_idx = i;
        else
            ID_labels(i) = 2;
            last_x_idx = i;
        end
        % second step of segment
        personID = impacts(i+1,4);
        if personID == 1 | personID == 2
            ID_labels(i+1) = 1;
            last_o_idx = i+1;
        else
            ID_labels(i+1) = 2;
            last_x_idx = i+1;
        end
        i = i + 2; % advance to 3rd step in segment
    else
        personID = impacts(i,4);
        if personID == 1 | personID == 2
            difference = impacts(i,1) - impacts(last_o_idx,1);
            ID_o = ID_o + difference/Fs_fsr;
            ID_o_count = ID_o_count + 1;
            last_o_idx = i;
            ID_labels(i) = 1;
        else
            difference = impacts(i,1) - impacts(last_x_idx,1);
            ID_x = ID_x + difference/Fs_fsr;
            ID_x_count = ID_x_count + 1;
            last_x_idx = i;
            ID_labels(i) = 2;
        end
    end
end

ID_o/ID_o_count
ID_x/ID_x_count

save(data_root_katie,'impactN','ID_o','ID_o_count','ID_x','ID_x_count',...
    'ID_labels','walk_edges','Fs_fsr','Fs_acc','-append')

%% step time ID cont. (doesn't work 3/22/22)
disp("new run")
estimateID = zeros(1,impactN);
confidence = zeros(1,impactN); % low is better!

step_o = 0.4675; % from above segment result
step_x = 0.6067;

% pre-load estimates with first two steps in every segment
for i = 1:impactN
    if walk_edges(i) == -1
        % just assign real values for now, figure this out later ##########
        estimateID(i) = ID_labels(i);
        estimateID(i+1) = ID_labels(i+1);
        i = i + 2;
    end
end
estimateID(21) = ID_labels(21); % bc first 2 of segment are both o
% SUM OF CONFIDENCE VALUES, MINIMIZE between all combinations?
% estimate
% 18, 38, 59, 79, 101, 122, 143, 164, 185, 206, 225
start_index = 39;
stop_index = 59;

% skip i values if going back (fix this when make this recursive)
redo = false;
redo_skipto_i = 0;

for i = start_index:stop_index %impactN % temporarily focus on individual walking segment
    if walk_edges(i) ~= 1 % if it's not the last step in the segment
        if walk_edges(i+1) ~= 1 % if it's not the second to last step
            if estimateID(i) ~= 0 % current node is assigned
                % rank next 3 in terms of step time
                currenti = impacts(i,1);
                if estimateID(i) == 1
                    curr_step = step_o;
                else
                    curr_step = step_x;
                end
                st1 = ( impacts(i+1,1) - currenti )/Fs_fsr;
                st2 = ( impacts(i+2,1) - currenti )/Fs_fsr;
                st3 = ( impacts(i+3,1) - currenti )/Fs_fsr;
                m1 = abs(st1-curr_step);
                m2 = abs(st2-curr_step);
                m3 = abs(st3-curr_step);
                m_arr = [m1,m2,m3];
                [score_val, score_i] = min(m_arr); % lowest score is highest confidence
                
                if estimateID(i+score_i) == 0 % estimated node is unassigned
                    estimateID(i+score_i) = estimateID(i);
                    confidence(i+score_i) = score_val;
                else % estimated node is already assigned
                    disp("here1")
                    i
                    score_i
                    m_arr
                    % use confidence to reassign if better; if not, assign
                    % 2nd highest score
                    if confidence(i+score_i) < score_val
                        % assign 2nd highest confidence
                        [sort_m_arr, sorted_i] = sort(m_arr);
                        if estimateID(i+sorted_i(2)) == 0 % if 2nd highest is not assigned
                            estimateID(i+sorted_i(2)) = estimateID(i);
                            confidence(i+sorted_i(2)) = sort_m_arr(2);
                        else
                            disp("here3")
                            i
                        end
                        if redo
                            redo = false;
                            i = redo_skipto_i;
                        end
                    elseif estimateID(i) ~= estimateID(i+score_i)
                        % reassign it, go back to the last time this was assigned and repeat
                        disp("here4")
                        confidence(i+score_i)
                        score_val
                        % find the last label
                        old_label = estimateID(i+score_i);
                        all_old_labels = find(estimateID(start_index:stop_index) == old_label);
                        last_assigned = all_old_labels(end-1) + start_index - 1;
                        redo = true;
                        redo_skipto_i = i+1;
                        
                        % reassign
                        estimateID(i+score_i) = estimateID(i);
                        confidence(i+score_i) = score_val;
                        
                        i = last_assigned
                    end
                end
            else % current node is unassigned
                disp("here2")
                i
            end
        elseif walk_edges(i+1) == 1 % it's the 2nd to last step
            if estimateID(i+1) == 0 % last step is not labeled
                estimateID(i+1) = estimateID(i);
            end
        elseif walk_edges(i+2) == 1 % it's the 3rd to last step
            % rank next 3 in terms of step time
                currenti = impacts(i,1);
                st1 = ( impacts(i+1,1) - currenti )/Fs_fsr;
                st2 = ( impacts(i+2,1) - currenti )/Fs_fsr;
                [m1,i1] = min([abs(st1-step_o), abs(st1-step_x)]);
                [m2,i2] = min([abs(st2-step_o), abs(st2-step_x)]);
                [score_val, score_i] = min([m1,m2]); % lowest score is highest confidence
                if estimateID(i+score_i) == 0 % estimated node is unassigned
                    estimateID(i+score_i) = estimateID(i);
                    confidence = score_val;
                else
                    disp("here3")
                    i
                end
        end
    end
end


figure;
real_o = find(ID_labels(start_index:stop_index) == 1);
real_i = find(ID_labels(start_index:stop_index) == 2);
title('Real ID labels')

plot(real_o,0,'bo')
hold on
plot(real_i,0,'bx')


est_o = find(estimateID(start_index:stop_index) == 1);
est_i = find(estimateID(start_index:stop_index) == 2);
figure;
plot(est_o,0,'bo')
hold on
plot(est_i,0,'bx')
title('Estimated ID labels')


%% look at all combinations, find one with min uncertainty 3/22/22

disp("new run")

% 18, 38, 59, 79, 101, 122, 143, 164, 185, 206, 225
segments = [0, 18, 38, 59, 79, 101, 122, 143, 164, 185, 206, 225];

for s = 2:length(segments)
    start_index = segments(s-1)+1;
    stop_index = segments(s);

    powerN = stop_index - start_index + 1;

    allcombinations = ones(1, powerN); % no x's

    for i = round(powerN/4):round(3*powerN/4) % know they at least have to be 25% of total footsteps
        curr_arr = ones(1,powerN);
        curr_arr(1:i) = 2;
        P = uniqueperms(curr_arr);
        allcombinations = [allcombinations; P];
    end

    allcombN = length(allcombinations);
    scores = zeros(allcombN,1);

    for i = 1:allcombN
        curr_arr = allcombinations(i,:);
        one_idx = find(curr_arr == 1);
        two_idx = find(curr_arr == 2);
        score = 0;
        for j = 2:length(one_idx)
            curr_i = one_idx(j) + start_index - 1;
            past_i = one_idx(j-1) + start_index - 1;
            curr_diff = abs(impacts(curr_i,1) - impacts(past_i,1))/Fs_fsr;
            score = score + abs(curr_diff - step_o);
        end
        for k = 2:length(two_idx)
            curr_i = two_idx(k) + start_index - 1;
            past_i = two_idx(k-1) + start_index - 1;
            curr_diff = abs(impacts(curr_i,1) - impacts(past_i,1))/Fs_fsr;
            score = score + abs(curr_diff - step_x);
        end
        scores(i) = score;
    end

    [~,idx] = min(scores);
    estimateID = allcombinations(idx,:);

    figure;
    real_o = find(ID_labels(start_index:stop_index) == 1);
    real_x = find(ID_labels(start_index:stop_index) == 2);
    plot(real_o,0,'bo')
    hold on
    plot(real_x,0,'bx')
    title(sprintf('Real ID labels %d segment'),s)


    est_o = find(estimateID == 1);
    est_x = find(estimateID == 2);
    figure;
    plot(est_o,0,'bo')
    hold on
    plot(est_x,0,'bx')
    title(sprintf('Estimated ID labels %d segment'),s)
    
    difference_results = ID_labels(start_index:stop_index) - estimateID;
    correct = find(abs(difference_results) == 0);
    accuracy = length(correct)/length(estimateID)
end

%% examining the errors from above 3/22/22

% 18, 38, 59, 79, 101, 122, 143, 164, 185, 206, 225
start_index = 207;
stop_index = 225;

powerN = stop_index - start_index + 1;

allcombinations = ones(1, powerN); % no x's

for i = round(powerN/4):round(3*powerN/4) % know they at least have to be 25% of total footsteps
    curr_arr = ones(1,powerN);
    curr_arr(1:i) = 2;
    P = uniqueperms(curr_arr);
    allcombinations = [allcombinations; P];
end

allcombN = length(allcombinations);
scores = zeros(allcombN,1);

for i = 1:allcombN
    curr_arr = allcombinations(i,:);
    one_idx = find(curr_arr == 1);
    two_idx = find(curr_arr == 2);
    score = 0;
    for j = 2:length(one_idx)
        curr_i = one_idx(j) + start_index - 1;
        past_i = one_idx(j-1) + start_index - 1;
        curr_diff = abs(impacts(curr_i,1) - impacts(past_i,1))/Fs_fsr;
        score = score + abs(curr_diff - step_o);
    end
    for k = 2:length(two_idx)
        curr_i = two_idx(k) + start_index - 1;
        past_i = two_idx(k-1) + start_index - 1;
        curr_diff = abs(impacts(curr_i,1) - impacts(past_i,1))/Fs_fsr;
        score = score + abs(curr_diff - step_x);
    end
    scores(i) = score;
end

[~,idx] = min(scores);
estimateID = allcombinations(idx,:);

figure;
real_o = find(ID_labels(start_index:stop_index) == 1);
real_x = find(ID_labels(start_index:stop_index) == 2);
plot(real_o,0,'bo')
hold on
plot(real_x,0,'bx')
title(sprintf('Real ID labels %d segment'),s)


est_o = find(estimateID == 1);
est_x = find(estimateID == 2);
figure;
plot(est_o,0,'bo')
hold on
plot(est_x,0,'bx')
title(sprintf('Estimated ID labels %d segment'),s)

difference_results = ID_labels(start_index:stop_index) - estimateID;
correct = find(abs(difference_results) == 0);
accuracy = length(correct)/length(estimateID)

step_times_abs = impacts(start_index:stop_index,1);
figure; 
plot(real_o, step_times_abs(real_o),'o')
hold on
plot(real_x, step_times_abs(real_x),'o')

figure;
plot(est_o, step_times_abs(est_o),'o')
hold on
plot(est_x, step_times_abs(est_x),'o')

%% finding impact starts from pcb data 3/23/22

clip_start = impacts(start_index,1)*Fs_pcb/Fs_fsr;
clip_end = impacts(stop_index,1)*Fs_pcb/Fs_fsr;
clip_impact_times = impacts(start_index:stop_index,1)*Fs_pcb/Fs_fsr;
clip_impact_times = (clip_impact_times - clip_impact_times(1))/Fs_pcb;
figure;
spectrogram(pcbData(clip_start:clip_end,1),256,250,[],Fs_pcb,'yaxis');
hold on
plot(clip_impact_times,0,'rx')

figure;
plot(pcbData(clip_start:clip_end,1))
hold on
plot(clip_impact_times*Fs_pcb,0,'rx')

% LPF ends up not helping
figure; lowpass(pcbData(clip_start:clip_end,1),100,Fs_pcb,'ImpulseResponse','iir','Steepness',0.8)

figure; bandpass(pcbData(clip_start:clip_end,1),[2000 3000], Fs_pcb,'ImpulseResponse','iir','Steepness',0.7)

%% cont.

% 18, 38, 59, 79, 101, 122, 143, 164, 185, 206, 225
start_time = impacts(1,1)/Fs_fsr - 0.25;
stop_time = impacts(18,1)/Fs_fsr + 0.25;

pcb_start_idx = round(start_time*Fs_pcb);
pcb_stop_idx = round(stop_time*Fs_pcb);

for i = 1:6
    figure;
%     plot(pcbData(:,i))
    [upper, lower] = envelope(pcbData(pcb_start_idx:pcb_stop_idx,i),300,'rms');
    plot(pcbTime(pcb_start_idx:pcb_stop_idx), upper)
    hold on
    plot(pcbTime(pcb_start_idx:pcb_stop_idx), lower)
    hold on
    plot(impacts(1:18,1)/Fs_fsr,0,'rx','MarkerSize',6)
    title(sprintf('Sensor %d',i))
end
%% svd doesn't work
svdM = pcbData(pcb_start_idx:pcb_stop_idx,:);
downsample_svdM = downsample(svdM,10);
[U,S,V] = svd(downsample_svdM.');

figure;
for i = 1:6
    subplot(6,1,i)
    plot(V(:,i))
end
title('V')

figure;
for i = 1:6
    subplot(6,1,i)
    plot(downsample_svdM(:,i))
end
title('raw data')

%% ica

[W,Zhat] = ica(transpose(downsample_svdM));
[Winv,Sinv] = jade(W,2);
[r,c] = size(Winv);

figure;
for i = 1:6
    subplot(6,1,i)
    plot(Zhat(i,:))
end

%% wiener filter doesn't work either

sample = zeros(1,length(downsample_svdM)); % make this 5 sec segment
increment = round(step_x*Fs_pcb/10); % downsampling rate
count = increment;
while count < length(downsample_svdM)
    sample(count) = 1;
    count = count + increment;
end
figure; plot(sample)
% wiener cont.
norm_raw = pcbData(113830:276028,1)./max(pcbData(113830:276028,1));
norm_sample = datas(260977:423175,1)./max(datas(260977:423175,1));
[xest,b,MSE] = wienerFilt(norm_raw, norm_sample ,200); % last is filter order

figure; plot(norm_raw)
title('raw')
figure; plot(xest)
title('xest')

%% scalogram
% figure;
[wt,f] = cwt(pcbData(113830:276028,1).*10,Fs_pcb);
t = linspace(0,162199/Fs_pcb,162199);
figure; contour(t,f,abs(wt))
set(gca, 'YScale', 'log')
%%

figure;
maxwt = zeros(1,length(f));
for i = 1:length(f)
    maxval = max(abs(wt(i,:)));
    maxwt(i) = maxval;
end
plot(maxwt)

valididx = find(f < 600 & f > 100);
sumwt = zeros(1,length(f));
for i = 1:length(f)
    sumval = max(abs(wt(i,valididx)));
    sumwt(i) = sumval;
end
figure;
plot(sumwt)
hold on
plot(impacttimes-pcbTime(113830),0,'rx','MarkerSize',8)

figure; contour(t,f(valididx),abs(wt(valididx,:))./max(abs(wt(valididx,:))))
hold on
plot(impacttimes-pcbTime(113830),150,'rx','MarkerSize',8)

figure;
plot(pcbTime(113830:276028)-pcbTime(113830) ,pcbData(113830:276028,1))
hold on
startimpact = 113830 * Fs_fsr / Fs_pcb;
lastimpact = 276028 * Fs_fsr / Fs_pcb;
impactidx = find(impacts(:,1) >= startimpact & impacts(:,1) <= lastimpact);
impacttimes = impacts(impactidx,1)/Fs_fsr;
plot(impacttimes-pcbTime(113830),0,'rx','MarkerSize',8)
%%
figure; helperPlotScalogram(pcbData(113830:276028,1),Fs_pcb)

%% 3/24/22 WORKS! for detecting footfalls

segments = [0, 18, 38, 59, 79, 101, 122, 143, 164, 185, 206, 225];
t = linspace(0,162199/Fs_pcb,162199);
figure;
plot(impacttimes-pcbTime(113830),150,'rx','MarkerSize',8)
hold on
for i = 1:6
    pcbclip = pcbData(113830:276028,i);
    [wt,f] = cwt(pcbclip, Fs_pcb);
    hold on
    contour(t,f(valididx),abs(wt(valididx,:)))
    hold on
    title(sprintf("Sensor %d",i))
%     set(gca, 'YScale', 'log')
end









