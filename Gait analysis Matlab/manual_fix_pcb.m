function [arrival_idx,peak_idx,peak_mag] = manual_fix_pcb(arrival_wrong,arrival_right,peak_idx_right,pcbTime,arrival_idx,peak_idx,peak_mag,filt_pcbD)
% manually fixes wrong arrival/peak times for pcb data
% cols are sensors 1-4
% NaN if only some of the sensors need to be changed
% 11/6/21

sizeM = size(arrival_wrong);
for i = 1:sizeM(1) % for # of rows in arrival_wrong:
    for j = 1:4 % change values per sensor
        if ~isnan(arrival_wrong(i,j))
            wrong_idx = findTindex(arrival_wrong(i,j),pcbTime);
            indeces = abs(arrival_idx(:,j) - wrong_idx);
            [~,change_idx] = min(indeces);
            arrival_idx(change_idx,j) = findTindex(arrival_right(i,j),pcbTime);
            peak_idx(change_idx,j) = findTindex(peak_idx_right(i,j),pcbTime);
            peak_mag(change_idx,j) = abs(filt_pcbD(peak_idx(change_idx,j),j));
        end
    end
end

figure;
for i=1:4
    subplot(4,1,i)
    plot(pcbTime,filt_pcbD(:,i))
    hold on
    plot(pcbTime(peak_idx(:,i)),peak_mag(:,i),'r.','MarkerSize',10)
    plot(pcbTime(arrival_idx(:,i)),0,'k.','MarkerSize',10)
end
    
end

