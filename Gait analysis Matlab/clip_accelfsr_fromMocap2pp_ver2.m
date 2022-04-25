function [accData, accTime, fsrData, fsrTime] = clip_accelfsr_fromMocap2pp_ver2(fsrData, accData, trigger_signal, Fs_fsr, Fs_acc, Fs_trigger)
% 4/18/22
% clips the fsr data to where trigger data peaks and dips to match mocap
% also resamples fsr and accel data
% different from w/o ver2 bc using new delsys software


idx = find(trigger_signal < -1);
starti = idx(1);
idx = find(trigger_signal(starti:end) > -1);
endi = idx(1) + starti - 1;

startt = starti/Fs_trigger;
endt = endi/Fs_trigger;

converted_starti_fsr = round(startt*Fs_fsr);
converted_endi_fsr = round(endt*Fs_fsr);

converted_starti_acc = round(startt*Fs_acc);
converted_endi_acc = round(endt*Fs_acc);

fsrTime = linspace(0,endt-startt,converted_endi_fsr-converted_starti_fsr+1);
accTime = linspace(0,endt-startt,converted_endi_acc-converted_starti_acc+1);
fsrData = fsrData(converted_starti_fsr:converted_endi_fsr,:);
accData = accData(converted_starti_acc:converted_endi_acc,:);

end