%% 6/17/20

% 1/13/21: go to Dropbox - analysis/data, "pwd" into command window to get
% this path
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\data\';
load([data_root_katie, 'data_01-10-2021_15-07-26'])
% load('data/data_01-10-2021_15-07-26') %'data/data_09-21-2020_15-29.mat'
loc_names = {'A1', 'A5', 'E1'};
% Ben test
Ben = 100;
Fs = 12800;
plotData(datas, times, loc_names, 'raw data ');
plotFreqSpectrum(datas, Fs, loc_names);
filt_datas = lpf_data(datas);
plotData(filt_datas, times, loc_names, 'filtered data ');
%test github
%% 6/19/20
[p_idx, p_val] = findPeaks(filt_datas, times, Fs, loc_names);

%% 6/22/20 manual cross correlation
% % for barefot/sneakers walking anywhere
% s_A1_segs = [21390, 24750, 28280,33160,36590,42140,46720,52250,57820,60800, ...
%     64010,67600,73610,77470,79990,82390,85920,91440,96460,101200,104500,110000];
% s_c = self_corr(s_A1, s_A1_segs);
% 
% b_A1_segs = [28690,33670,36580,41600,43760,50670,59150,65120,69680,76040,...
%     79190,84580,88770,93460];
% b_c = self_corr(b_A1, b_A1_segs);
% 
% sb_c = comp_corr(b_A1, b_A1_segs, s_A1, s_A1_segs);

% for barefoot/sneakers walking D3 to B3
% load('data\data_06-23-2020_10-50')
% filt_b = lpf_data(datas);
% b_A5 = filt_b(:,3);
% load('data\data_06-23-2020_10-52')
% filt_s = lpf_data(datas);
% s_A5 = filt_s(:,3);
% 
% s_A5_segs = [43210,61220,77190,94690,113700,130900,152800,164800,184900,...
%     201500,218700,236400,263000,274800,297400,315800];
% b_A5_segs = [42430,55450,70840,88110,102900,118100,133500,147200,163600,...
%     178800,193700,206900];
% 
% s_c = self_corr(s_A5, s_A5_segs);
% b_c = self_corr(b_A5, b_A5_segs);
% sb_c = comp_corr(b_A5, b_A5_segs, s_A5, s_A5_segs);

% for barefoot/sneakers walking B5, B3, B1 & back
load('data/data_06-23-2020_15-24')
filt_b = lpf_data(datas);
b_A5 = filt_b(:,3);
load('data/data_06-23-2020_15-27')
filt_s = lpf_data(datas);
s_A5 = filt_s(:,3);

s_A5_segs = [98460,102100,266800,270700];
b_A5_segs = [109700,111300,298100,299400];

s_c = self_corr(s_A5, s_A5_segs);
b_c = self_corr(b_A5, b_A5_segs);
sb_c = comp_corr(b_A5, b_A5_segs, s_A5, s_A5_segs);

%% 6/29/20 localization B5,B3,B1
% M1: onset of impacts at B5, B3, B1, B1, B3, B5 (first large disturbance)
% first row is A1, second E1, then A5
M1_barefoot = [110500, 113200, 118700, 208700, 213100, 217600;
    110400, 113200, 118900, 209100, 213100, 216700;
    110400, 113200, 118900, 209100, 213100, 217700];
% M2: same as M1 but second iteration
M2_barefoot = [298300, 300800, 306100, 388900, 395400, 401100;
    298200, 300700, 306200, 389100, 395400, 400600;
    298300, 300700, 306200, 389200, 395400, 401100];
% M1_sneakers = [938700, 99000, 102400, 177200, 180100, 185400;
%     938500, 99000, 102400

% corners_middle_twice_07-02-2020_13-52
% these are the onset times of each impact, E1,E5,C3,A5,A1, twice each
M = [11.1871,15.8796,22.0016,27.0015,35.6155,40.0005,49.414,54.6445,59.7037,64.0013;...
    11.1859,15.8797,22.0012,27.0023,35.6151,40.0006,49.412,54.642,59.7065,64.0026;...
    11.1822,15.8763,22.0026,27.0022,35.6151,40.0013,49.415,54.6446,59.7066,64.0015];
Amp = [0.3148,0.2848,0.1651,0.05117,0.4521,0.08451,0.4729,0.5178,3.079,0.09194;...
    0.1654,0.1715,0.5543,0.06761,1.12,0.1161,6.492,5.757,0.2753,0.08903;...
    4.925,5.318,0.1039,0.01785,0.2298,0.03659,0.1157,0.1318,0.1859,0.03102];
[mag_diff, diff] = triangulate_timing(transpose(M),transpose(Amp), loc_names);


%% 7/3/20 TDOA on clean envelope data corners middle
load('data/corners_middle_twice_07-02-2020_13-52_07-02-2020_13-52.mat')
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);
% corners_middle_twice_07-02-2020_13-52
% these are the onset times of each impact, E1,E5,C3,A5,A1, twice each
M = [95443,135498,187742,230401,303902,341329,421649,466274,509454,546122;...
    95449,135502,187732,230411,303906,341329,421632,466262,509476,546135;...
    95409,135474,187748,230408,303907,341333,421659,466284,509478,546125];
Amp = [0.3266,0.2983,0.1651,0.05117,0.4521,0.08451,0.4729,0.5178,3.079,0.1145;...
    0.2131,0.2,0.5991,0.07248,1.12,0.1149,6.492,5.757,0.3349,0.1026;...
    4.925,5.318,0.1041,0.02407,0.2298,0.04326,0.1622,0.1774,0.1859,0.03465];
[mag_diff, diff] = triangulate_timing(transpose(M), transpose(Amp), loc_names);

%% 7/9/20 TDOA on clean envelope data corners middle
load('data/data_07-09-2020_13-11.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);
plotData(clean_data, times, loc_names, 'clean filtered data ');
% these are the onset times of each impact, E1,E3,E5,C5,C1,C3,A5,A1, twice each
M = [211715,286300,374854,450232,544553,624565,743568,813543,931214,1017843,1266236,1367050,1479860,1553015,1691914,1779470;...
    211721,286305,374864,450239,544552,624572,743561,813533,931224,1017853,1266237,1367050,1479844,1552999,1691937,1779492;...
    211685,286269,374838,450218,544544,624563,743565,813536,931214,1017845,1266232,1367048,1479870,1553027,1691939,1779495];
Amp = [0.3040485,0.330426,0.2429059,0.2625341,0.4017302,0.3744977,0.2347931,0.2469461,0.8321206,0.5943623,0.4536783,0.4899976,0.64819,0.8751741,3.870439,6.335793;...
    0.1664618,0.1466173,0.6883936,0.7546019,0.8254613,0.7441817,1.361391,1.687344,0.3799721,0.4481337,0.9788601,0.7947675,5.428729,5.388199,0.3388803,0.436588;...
    1.651993,4.058874,0.938072,0.9614384,0.6631244,0.5237078,0.1539047,0.1504284,0.9799187,1.103061,0.195569,0.2299322,0.1688069,0.1692021,0.1986951,0.2090298];
[mag_diff, diff] = triangulate_timing(transpose(M), transpose(Amp), loc_names);

%% 7/10/20 TDOA on clean envelope data corners middle
load('data/data_07-10-2020_14-40.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);
plotData(clean_data, times, loc_names, 'clean filtered data ');

M = [182522,233834,289997,339071,397803,441102,486337,527406,572122,612252,654057,694406,739934,779640,818904,854792;...
    182548,233840,290013,339092,397820,441052,486259,527421,572081,612257,654033,694387,739916,779550,818924,854803;...
    182510,233834,289938,339056,397724,441025,486213,527398,572090,612269,654068,694410,739939,779650,818947,854808];
Amp = [0.4703078,0.8401363,0.230835,0.2388231,0.1631428,0.1646059,0.3424525,0.1492062,0.1152247,0.097679,0.2258736,0.3742097,0.212593,0.5311266,0.8254342,3.244074;...
    0.2311123,0.4853187,0.1714903,0.1790185,0.4945618,0.304286,0.7609695,0.4305096,0.6040992,0.464942,1.828834,1.847711,2.129028,1.986123,0.2319196,0.2336708;...
    0.5589042,0.8876641,0.7454762,1.604415,0.8001235,0.8498281,0.3718151,0.1792815,0.1155433,0.0529052,0.09858681,0.106749,0.09146868,0.1395859,0.1760163,0.135868];
M_peaks = [182875,233938,290182,339209,398272,441445,486583,527819,572243,612538,654410,694838,740367,779767,819115,854974;...
    182813,234141,290247,339355,398134,441414,486676,527912,572274,612546,654228,694664,740416,779911,819161,855153;...
    182778,233998,289987,339177,398013,441242,486468,527788,572255,612637,654409,694840,740464,779932,819256,855053];
[mag_diff, diff] = triangulate_timing(transpose(M), transpose(Amp), loc_names);
% [mag_diff2, diff2] = triangulate_timing(transpose(M_peaks), transpose(Amp), loc_names);

%% 7/24/20 algorithm for TDOA on ball bounce
clear all
load('data/data_10-02-2020_19-03.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

impactN = 16; % variable input value

[onset_idx, peak_idx, peak_val] = TDOA(clean_data,impactN,Fs,loc_names);
[mag_diff, diff] = triangulate_timing(onset_idx, peak_val, loc_names,Fs);


%% 7/28/20 kmeans localization

load('data/data_07-27-2020_14-26.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
% filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

impactN = 15; % variable input value

[onset_idx, peak_idx, peak_val] = TDOA(clean_data,impactN,Fs,loc_names);
[mag_diff, diff] = triangulate_timing(onset_idx, peak_val, loc_names,Fs);

dataset = addData([],diff,mag_diff);

load('data/data_07-27-2020_14-32.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

impactN = 15; % variable input value

[onset_idx, peak_idx, peak_val] = TDOA(clean_data,impactN,Fs,loc_names);
[mag_diff, diff] = triangulate_timing(onset_idx, peak_val, loc_names,Fs);

dataset = addData(dataset,diff,mag_diff);

idx = kmeans(dataset,6);

%% writing it to csv file for hw2 2.986 9/16/20
% data_07-27-2020_14-26.mat & data_07-27-2020_14-22.mat
d_names = ["Time delay A1", "Time delay A5", "Time delay E1", "Voltage offset A1", "Voltage offset A5", "Voltage offset E1", "Labels"];
labels = [1,1,1,1,1,3,3,3,3,3,5,5,5,5,5,1,1,1,1,1,3,3,3,3,3,5,5,5,5,5];
datasetL = [dataset,transpose(labels)];
datasetLN = [d_names;datasetL];
writematrix(datasetLN, 'dataset.csv');

%% gait imbalance recognition 8/3/20

load('data/data_09-25-2020_16-01.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

impactN = 13; % variable input value

[onset_idx, peak_idx, peak_val] = TDOA(clean_data,impactN,Fs,loc_names);
[mag_diff, diff] = triangulate_timing(onset_idx, peak_val, loc_names,Fs);
[A,B] = gaitAnalysis(onset_idx,Fs,impactN);

%% 10/1/20 testing new peak algorithm (TDOA2)
clear all
load('data/data_10-19-2020_15-49.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

impactN = 6; % ADJUST number of expected impacts in data

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
[mag_diff, diff] = triangulate_timing(transpose(onset_idx), transpose(peak_val), loc_names,Fs);
