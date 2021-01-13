load('data\data_07-09-2020_13-11.mat')
loc_names = {'A1', 'A5', 'E1'};

Fs = 12800;
plotData(datas, times, loc_names, 'raw data ');
plotFreqSpectrum(datas, Fs, loc_names);
filt_datas = lpf_data(datas);
plotData(filt_datas, times, loc_names, 'filtered data ');
clean_data = clean_envelope(filt_datas,Fs);
plotData(clean_data, times, loc_names, 'clean filtered data ');

impactN = 16; % variable input value
[onset_idx, peak_idx, peak_val] = TDOA(clean_data,impactN,Fs,loc_names);
[mag_diff, diff] = triangulate_timing(onset_idx, peak_val, loc_names);