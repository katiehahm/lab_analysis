load('data/data_10-02-2020_19-03.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

impactN = 16; % variable input value

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
[peakdistdifftable] = stridetime_anaylsisfunction(peak_idx,impactN);


