function [pcbData, pcbTime] = clip_pcb_fromMocap(datas, times)
% 8/30/21
% clips the pcb data where 6th channel peaks and dips to match mocap
idx = find(datas(:,8) > 1);
starti = idx(1);
idx = find(datas(starti:end,8) < 1);
endi = idx(1) + starti - 1;
pcbData = datas(starti:endi,1:4);
pcbTime = times(1:(endi-starti+1));

end

