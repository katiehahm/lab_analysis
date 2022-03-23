function [pcbData, pcbTime] = clip_pcb_fromMocap(datas, times, sensorN)
% 8/30/21
% clips the pcb data where 6th channel peaks and dips to match mocap
idx = find(datas(:,end) > 1);
starti = idx(1);
idx = find(datas(starti:end,end) < 1);
endi = idx(1) + starti - 1;
pcbData = datas(starti:endi,1:sensorN);
pcbTime = times(1:(endi-starti+1));

end

