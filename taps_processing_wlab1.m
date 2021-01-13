%% 3/3/20 post lab1
clear all
close all
load('footsteps_013020')
load('multitaps')

Fs = 192e+3;
d = balls;
% d = multitap_B4_E5;
% d = multitap_all_A4F6;

% adjust data so it is visually corresponding to sensor quadrants
dRT = d(:,RT);
dLT = d(:,LT);
dRB = d(:,RB);
dLB = d(:,LB);
d(:,1) = dLT;
d(:,2) = dRT;
d(:,3) = dLB;
d(:,4) = dRB;
loc_names = {'left top', 'right top', 'left bottom', 'right bottom'};
n = length(dRB);

figure; % raw data
hold on
for i = 1:4
    subplot(2,2,i)
    title_str = 'raw data ';
    plot(d(:,i))
    title([title_str, loc_names(i)])
end

figure; % frequency spectrum
hold on
for i = 1:4
    subplot(2,2,i)
    title_str = 'frequency spectrum ';
    [X, F] = pwelch(d(:,i),[],[],[],Fs);
    plot(F, pow2db(X))
    title([title_str, loc_names(i)])
end

%% filter
filt = zeros(n, 4);
for i = 1:4
%     b = fir1(48, 20000/Fs, 'low');
%     filt(:,i) = filter(b, 1, d(:,i));
    filt(:,i) = lowpass(d(:,i), 10000, Fs);
end

figure; % filtered signal in time domain
hold on
for i = 1:4
    subplot(2,2,i)
    title_str = 'filtered data ';
    plot(filt(:,i))
    title([title_str, loc_names(i)])
end

%% find peaks
[p_idx, p_val] = findPeaks(filt, Fs);

%% find onset of first disturbance
% [disturb, disturb_val] = onset(p_idx, filt);
% figure;
% hold on
% for i = 1:4
%     subplot(2,2,i)
%     title_str = 'filtered data with disturbance onsets ';
%     plot(filt(:,i))
%     hold on
%     plot(disturb(:,i), disturb_val(:,i), 'ro')
%     title([title_str, loc_names(i)])
% end
%% perform triangulation on peaks
[mag_diff, diff] = triangulate_timing(p_idx, p_val);

figure;
bar(diff)
legend('right top', 'right bottom', 'left top', 'left bottom')
title('Triangulation: time delays of each sensor from peaks')
xlabel('Hit')
ylabel('samples (*1000/Fs = ms)')

figure;
bar(mag_diff)
legend('right top', 'right bottom', 'left top', 'left bottom')
title('Triangulation: magnitude differences of each sensor at peak')
xlabel('Hit')
ylabel('Volts (+50 dB gain)')

%% perform triangulation on onset

% [mag_diff, diff] = triangulate_timing(disturb, disturb_val);
% 
% figure;
% bar(diff)
% legend('right top', 'right bottom', 'left top', 'left bottom')
% title('Triangulation: time delays of each sensor from onset')
% xlabel('Hit')
% ylabel('samples (*1000/Fs = ms)')
% 
% figure;
% bar(mag_diff)
% legend('right top', 'right bottom', 'left top', 'left bottom')
% title('Triangulation: magnitude differences of each sensor at onset')
% xlabel('Hit')
% ylabel('Volts (+50 dB gain)')



%% Guess from "similarity" from bank

load('peaksBank')
nHits = length(diff);
name_labels = [' A4'; ' A5'; ' A6'; ' B4'; ' B5'; ' B6'; ' C4'; ' C5'; ' C6'; ' D4'; ' D5'; ' D6'; ' E4'; ' E5'; ' E6'; ' F4'; ' F5'; ' F6'];
labels = '';
for i = 1:nHits
    diff_array = diff(i,:);
    pk_array = peaksBank(1,:);
    comp_array = pk_array(1:4);
    similarity = sum(abs(diff_array-comp_array));
    label = 1;
    for j = 2:length(peaksBank)
        pk_array = peaksBank(j,:);
        comp_array = pk_array(1:4); % take out label number
        new_diff = sum(abs(diff_array-comp_array));
        if new_diff < similarity
            similarity = new_diff;
            label = pk_array(end);
        end
    end
    labels = append(labels, name_labels(label,:));
end
labels
% find differences
% compare differences to mean of each square
% return location in strings


































