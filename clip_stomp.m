function [clip_stomp_data,clip_stomp_time,clip_fsrD,clip_fsrT] = clip_stomp(clipped_time,clipped_data,fsrT,fsrD)
% 5/27/21
% clip both fsr and pcb data to make the big stomp t=0

[maxv,~] = max(clipped_data(1:50000,2));
datalen = length(clipped_data);
buffer = 0.5;
for i = 1:datalen
   if clipped_data(i, 2) > maxv + buffer % always stomp next to S2
       break
   end
end

window = clipped_data(i-3000:i+10000,2); % approx. before and after impact
arrival_idx = aic_pick(window, 'to_peak')+i-3000-1;
clipped_time(arrival_idx) % this is the time pcb thinks impact happened

fsri = findTindex(clipped_time(arrival_idx),fsrT);
Rheeli = 3;
window2 = fsrD(Rheeli, fsri- 200:fsri+400);
arrival_idx2 = aic_pick(window2, 'to_peak')+fsri-200-1;
fsrT(arrival_idx2) % this value should be similar to above

clip_stomp_data = clipped_data(arrival_idx:end,:);
clip_stomp_time = linspace(0,clipped_time(end)-clipped_time(arrival_idx), length(clipped_time)-arrival_idx+1);

clip_fsrD = fsrD(:,arrival_idx2:end);
clip_fsrT = linspace(0, fsrT(end)-fsrT(arrival_idx2), length(fsrT)-arrival_idx2+1);

end

