function [] = save_processed_data(filename, n, sn, m)
% 1/24/21
% to run on all data files to save features into the file
% to be then processed (to not duplicate feature finding)

% clean_data = clean_envelope(filt_data,Fs); % took abs val out 1/24/21

% changeable parameters
impactN = n;
impactWidth = 10000; % number of indeces an impact lasts in data
mass = m; % mass of object if dropped

arrival_idx = zeros(n,sn);
peak_idx = zeros(n,sn);
peak_mag = zeros(n,sn);
prominence = zeros(n,sn);
echo = zeros(n,sn);
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
        
        prominence(i,j) = peak_mag(i,j) - min(window);
        echo(i,j) = sum(abs(window)); % Can get rid of abs but this was how I originally did 
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

<<<<<<< Updated upstream
save([data_root_katie, filename], 'raw_data','loc_names','Fs','filt_data','impactN','arrival_idx','peak_idx','peak_mag', 'prominence', 'echo', 'mass', '-append')
=======
% save([data_root_katie, filename], 'raw_data','loc_names','Fs','filt_data','impactN','arrival_idx','peak_idx','peak_mag','mass', '-append')
>>>>>>> Stashed changes
   
end

