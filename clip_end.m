function [pcbD,pcbT,fsrT,fsrD,mocapT,mocapR,mocapL] = clip_end(clip_stomp_data,clip_stomp_time,clip_fsrT,clip_fsrD,mocapT,mocapR,mocapL)
% 5/27/21
% clip the end of the data on all 3 sensors 
% this is mostly walking back to the computer
% this clips after the floor sensor S2 don't hear anything anymore.

[maxv,~] = max(clip_stomp_data(end-100:end,:));

datalen = length(clip_stomp_data);
for i = 1:datalen
    if any(clip_stomp_data(datalen+1-i,:) > maxv+0.05)
        break
    end
end

buffer = 100;
lastindex = datalen+1-i + buffer;

pcbD = clip_stomp_data(1:lastindex,:);
pcbT = clip_stomp_time(1:lastindex);

lastindexfsr = findTindex(pcbT(end),clip_fsrT);
fsrT = clip_fsrT(1:lastindexfsr);
fsrD = clip_fsrD(:,1:lastindexfsr);

lastindexmc = findTindex(pcbT(end),mocapT);
mocapL = mocapL(1:lastindexmc,:);
mocapR = mocapR(1:lastindexmc,:);
mocapT = mocapT(1:lastindexmc);

end

