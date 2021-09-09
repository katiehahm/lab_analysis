function [fsrData, fsrTime] = clip_fsr_fromMocap(Data, Time)
% 8/30/21
% clips the fsr data where 9th channel peaks and dips to match mocap
% also resamples data: FSR is 519 Hz, Trigger is 2222 Hz
idx = find(Data(9,:) < -1);
starti = idx(1);
idx = find(Data(9,starti:end) > -1);
endi = idx(1) + starti - 1;

startt = Time(9,starti);
endt = Time(9,endi);

converted_starti = round(startt*519);
converted_endi = round(endt*519);

fsrTime = linspace(0, Time(9,end), Time(9,end)*519);
Data = Data(1:8,converted_starti:converted_endi);
fsrData = Data';

fsrTime = fsrTime(1:(converted_endi-converted_starti+1));


end