% creating dataset for quiz1 2s986
% 10/8/20


keySet = [1,2,3,4,5,6];
valueSet = ['BQ', 'CQ', 'DQ', 'EQ', 'ER', 'ES'];
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
%%
load('data/data_09-25-2020_15-27.mat')
filt_datas = lpf_data(datas);

filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

impactN = 10;

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
[mag_diff, diff] = triangulate_timing(transpose(onset_idx), transpose(peak_val), loc_names,Fs);
labels = [1,1,1,1,1,2,2,2,2,2];
datatable = table(mag_diff, diff, transpose(labels));

load('data/data_09-25-2020_15-28.mat')
filt_datas = lpf_data(datas);
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

impactN = 10;

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
[mag_diff, diff] = triangulate_timing(transpose(onset_idx), transpose(peak_val), loc_names,Fs);
labels = [3,3,3,3,3,4,4,4,4,4];
dataT = table(mag_diff, diff, transpose(labels));
datatable = [datatable; dataT];

load('data/data_09-25-2020_15-29.mat')
filt_datas = lpf_data(datas);
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

impactN = 10;

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
[mag_diff, diff] = triangulate_timing(transpose(onset_idx), transpose(peak_val), loc_names,Fs);
labels = [5,5,5,5,5,6,6,6,6,6];
dataT = table(mag_diff, diff, transpose(labels));
datatable = [datatable; dataT];

load('data/data_09-25-2020_15-30.mat')
filt_datas = lpf_data(datas);
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

impactN = 10;

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
[mag_diff, diff] = triangulate_timing(transpose(onset_idx), transpose(peak_val), loc_names,Fs);
labels = [1,1,1,1,1,2,2,2,2,2];
dataT = table(mag_diff, diff, transpose(labels));
datatable = [datatable; dataT];

load('data/data_09-25-2020_15-33.mat')
filt_datas = lpf_data(datas);
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

impactN = 10;

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
[mag_diff, diff] = triangulate_timing(transpose(onset_idx), transpose(peak_val), loc_names,Fs);
labels = [3,3,3,3,3,4,4,4,4,4];
dataT = table(mag_diff, diff, transpose(labels));
datatable = [datatable; dataT];

load('data/data_09-25-2020_15-31.mat')
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

impactN = 10;

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
[mag_diff, diff] = triangulate_timing(transpose(onset_idx), transpose(peak_val), loc_names,Fs);
labels = [5,5,5,5,5,6,6,6,6,6];
dataT = table(mag_diff, diff, transpose(labels));
datatable = [datatable; dataT];
%%
load('data/data_09-18-2020_11-22.mat')
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

impactN = 10;

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
[mag_diff, diff] = triangulate_timing(transpose(onset_idx), transpose(peak_val), loc_names,Fs);
labels = [2,2,2,2,2,2,2,2,2,2];
dataT = table(mag_diff, diff, transpose(labels));
datatable = [datatable; dataT];

load('data/data_09-18-2020_11-23.mat')
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

impactN = 10;

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
[mag_diff, diff] = triangulate_timing(transpose(onset_idx), transpose(peak_val), loc_names,Fs);
labels = [2,2,2,2,2,2,2,2,2,2];
dataT = table(mag_diff, diff, transpose(labels));
datatable = [datatable; dataT];