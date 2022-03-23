function [accData, accTime, fsrData, fsrTime] = clip_accelfsr_fromMocap2pp(Data, Time, Fs, TriggerRow)
% 11/22/21
% clips the fsr data where 9th channel peaks and dips to match mocap
% also resamples data: FSR is 296.3 Hz, Trigger is 2222 Hz
% also resamples accel data: 148.1 Hz
% Fs is the fsr data sampling rate. Trigger is 2222.2 Hz. 

Fs_heel = Fs(1);
Fs_acc = Fs(2);

idx = find(Data(TriggerRow,:) < -1);
starti = idx(1);
idx = find(Data(TriggerRow,starti:end) > -1);
endi = idx(1) + starti - 1;

startt = Time(TriggerRow,starti);
endt = Time(TriggerRow,endi);

% delete trigger row from data so that the fsr/accel row numbers are
% consistent
Data(TriggerRow,:) = [];

converted_starti_heel = round(startt*Fs_heel);
converted_endi_heel = round(endt*Fs_heel);

converted_starti_acc = round(startt*Fs_acc);
converted_endi_acc = round(endt*Fs_acc);

maxT = max(Time(TriggerRow,:)); % because sometimes end of Time is 0
fsrTime = linspace(0, maxT, maxT*Fs_heel);
fsrTime = fsrTime(1:(converted_endi_heel-converted_starti_heel+1));
fsrData = Data(1,converted_starti_heel:converted_endi_heel);
fsrData = [fsrData; Data(5,converted_starti_heel:converted_endi_heel)];
fsrData = [fsrData; Data(9,converted_starti_heel:converted_endi_heel)];
fsrData = [fsrData; Data(13,converted_starti_heel:converted_endi_heel)];
fsrData = fsrData';

accTime = linspace(0, maxT, maxT*Fs_acc);
accTime = accTime(1:(converted_endi_acc-converted_starti_acc+1));
accData = Data(2:4,converted_starti_acc:converted_endi_acc);
accData = [accData; Data(6:8,converted_starti_acc:converted_endi_acc)];
accData = [accData; Data(10:12,converted_starti_acc:converted_endi_acc)];
accData = [accData; Data(14:16,converted_starti_acc:converted_endi_acc)];
accData = accData';

end