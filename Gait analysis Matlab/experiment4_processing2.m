close all
clear all
processedfilepath = 'C:\Users\katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Jenny 1\ProcessedData\both_limp2.mat';
load(processedfilepath)
% make training set 6/10/22
freq_lower = 100; % Hz; frequency limits for cwt
freq_higher = 350;
arrive_idx_all = zeros(impactN,sensorN);
last_idx_all = zeros(impactN,sensorN);
width = zeros(impactN,sensorN);
peak_mag = zeros(impactN,sensorN);
energy = zeros(impactN,sensorN);
cwt_energy = zeros(impactN,sensorN);
cwt_peak = zeros(impactN,sensorN);
overlapping = zeros(impactN,1);
overlap_thresh = 1000; % # indexes that's considered an overlapping impact
window_width = 0.35; % in (s) of how long to set initial lastidx of window

Y = impacts(:,4);
Y(Y==2) = 1;
Y(Y==3) = 2;
Y(Y==4) = 2;

for seg = 2:length(segments)
    start_impact_idx = segments(seg-1);
    last_impact_idx = segments(seg)-1;
    start_time_fsr = impacts(start_impact_idx:last_impact_idx,1)./Fs_fsr;
    pcb_start_idx = round((start_time_fsr(1)-0.5)*Fs_pcb);
    pcb_last_idx = min(round((start_time_fsr(end)+0.5)*Fs_pcb),length(wien_pcbD(:,1)));
    for s = 1:sensorN
        seg_pcb = wien_pcbD(pcb_start_idx:pcb_last_idx,s);
        [wt,f] = cwt(seg_pcb,Fs_pcb);
        valid_f_idx = find(freq_higher & f > freq_lower);
        cwt_mag = abs(wt(valid_f_idx,:));
        sum_cwt = sum(cwt_mag,1);
        sum_smooth_cwt = movmean(sum_cwt, 1000);
        figure; 
        plot(sum_smooth_cwt)
        hold on
        for i = 1:length(start_time_fsr)
            global_idx = start_impact_idx + i - 1;
            fsr_start_idx = round(start_time_fsr(i)*Fs_pcb);
            starti = fsr_start_idx - 1000;
            if i ~= length(start_time_fsr) % not last element
                if start_time_fsr(i+1)*Fs_pcb - start_time_fsr(i)*Fs_pcb > overlap_thresh
                    lasti = min([round(start_time_fsr(i+1)*Fs_pcb), starti + round(window_width*Fs_pcb)]);
                else
                    lasti = starti + round(window_width*Fs_pcb);
                end
            else
                lasti = min([pcb_last_idx,starti + round(window_width*Fs_pcb)]);
            end
            cwt_window = sum_smooth_cwt(starti-pcb_start_idx:lasti-pcb_start_idx);
            % find arrival idx with aic pick
            aic_idx = aic_pick(cwt_window, 'to_peak');
            arrive_idx = aic_idx + starti;
            if arrive_idx < fsr_start_idx
                arrive_idx = fsr_start_idx;
            end
            % last_idx is the 1st idx that goes under thresh after peak
            cwt_window = sum_smooth_cwt(arrive_idx:lasti-pcb_start_idx);
            [~,maxidx] = max(cwt_window);
            under_thresh_idx = find(cwt_window(maxidx:end) < 0.00065);
            if isempty(under_thresh_idx)
                last_idx = lasti;
            else
                last_idx = under_thresh_idx(1) + aic_idx - 1 + maxidx - 1 + starti;
            end
            % store vars
            arrive_idx_all(global_idx,s) = arrive_idx;
            last_idx_all(global_idx,s) = last_idx;
            width(global_idx,s) = last_idx - arrive_idx;
            cwt_window = sum_smooth_cwt(arrive_idx-pcb_start_idx:last_idx-pcb_start_idx);
            cwt_peak(global_idx,s) = max(cwt_window);
            cwt_energy(global_idx,s) = sum(abs(cwt_window).^2);
            pcb_window = wien_pcbD(arrive_idx:last_idx,s);
            peak_mag(global_idx,s) = max(pcb_window);
            energy(global_idx,s) = sum(abs(pcb_window).^2);
            % overlap check
            if i > 1
                if round(start_time_fsr(i)*Fs_pcb) - round(start_time_fsr(i-1)*Fs_pcb) < 1000
                    overlapping(global_idx) = 1;
                    overlapping(global_idx-1) = 1;
                end
            end
            % sanity check plot
            plot(arrive_idx-pcb_start_idx,0,'rx','MarkerSize',8)
            plot(last_idx-pcb_start_idx,0,'gx','MarkerSize',12)
        end
        
    end
end

X = [width,cwt_peak,cwt_energy,peak_mag,energy];
overlap_idx = find(overlapping == 1);
Y(overlap_idx) = [];
X(overlap_idx,:) = [];

X_train = X;
Y_train = Y;

% make overall impacts variable
old_impacts = impacts;
person_labeling = old_impacts(:,4);
person_labeling(person_labeling == 1) = 11;
person_labeling(person_labeling == 2) = 12;
person_labeling(person_labeling == 3) = 21;
person_labeling(person_labeling == 4) = 22;

impacts = [impacts(:,1)./Fs_fsr,person_labeling, coordinates(:,1),coordinates(:,3),acc_pks(:,2),overlapping];

%% save vars
save(processedfilepath,'arrive_idx_all','last_idx_all','width','peak_mag','energy',...
    'cwt_energy','cwt_peak','overlapping','X_train','Y_train',...
    'old_impacts','impacts','-append')












