function [fsrData, fsrTime] = clip_fsr_fromMocap(Data, Time, Fs)
% 8/30/21
% clips the fsr data where 9th channel peaks and dips to match mocap
% also resamples data: FSR is 518.5 Hz, Trigger is 2222 Hz
% Fs is the fsr data sampling rate. Trigger is 2222 Hz. 
lastI = length(Data(:,1));
idx = find(Data(lastI,:) < -1);
starti = idx(1);
idx = find(Data(lastI,starti:end) > -1);
endi = idx(1) + starti - 1;

startt = Time(lastI,starti);
endt = Time(lastI,endi);

converted_starti = round(startt*Fs);
converted_endi = round(endt*Fs);

maxT = max(Time(lastI,:)); % because sometimes end of Time is 0
fsrTime = linspace(0, maxT, maxT*Fs);
Data = Data(1:(lastI-1),converted_starti:converted_endi);
fsrData = Data';

fsrTime = fsrTime(1:(converted_endi-converted_starti+1));


end