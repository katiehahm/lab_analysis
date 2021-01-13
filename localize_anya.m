clear all
load('data/data_12-10-2020_12-16-21.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1', 'E5'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope2(filt_datas,Fs);

impactN = 15; % ADJUST number of expected impacts in data

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
avgs1=[];
avgs2=[];
avgs3=[];
avgs4=[];
for i=1:5:15
    avg1= mean(peak_val(1,i:i+4));
    avgs1=[avgs1 avg1];
    avg2= mean(peak_val(2,i:i+4));
    avgs2=[avgs2 avg2];
    avg3= mean(peak_val(3,i:i+4));
    avgs3=[avgs3 avg3];
    avg4= mean(peak_val(4,i:i+4));
    avgs4=[avgs4 avg4];
end

load('data/data_12-10-2020_12-18-11.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1', 'E5'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope2(filt_datas,Fs);

impactN = 15; % ADJUST number of expected impacts in data

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);

for i=1:5:15
    avg1= mean(peak_val(1,i:i+4));
    avgs1=[avgs1 avg1];
    avg2= mean(peak_val(2,i:i+4));
    avgs2=[avgs2 avg2];
    avg3= mean(peak_val(3,i:i+4));
    avgs3=[avgs3 avg3];
    avg4= mean(peak_val(4,i:i+4));
    avgs4=[avgs4 avg4];
end

load('data/data_12-10-2020_12-19-54.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1', 'E5'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope2(filt_datas,Fs);

impactN = 15; % ADJUST number of expected impacts in data

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);

for i=1:5:15
    avg1= mean(peak_val(1,i:i+4));
    avgs1=[avgs1 avg1];
    avg2= mean(peak_val(2,i:i+4));
    avgs2=[avgs2 avg2];
    avg3= mean(peak_val(3,i:i+4));
    avgs3=[avgs3 avg3];
    avg4= mean(peak_val(4,i:i+4));
    avgs4=[avgs4 avg4];
end

load('data/data_12-10-2020_12-21-17.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1', 'E5'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope2(filt_datas,Fs);

impactN = 10; % ADJUST number of expected impacts in data

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);

for i=1:5:10
    avg1= mean(peak_val(1,i:i+4));
    avgs1=[avgs1 avg1];
    avg2= mean(peak_val(2,i:i+4));
    avgs2=[avgs2 avg2];
    avg3= mean(peak_val(3,i:i+4));
    avgs3=[avgs3 avg3];
    avg4= mean(peak_val(4,i:i+4));
    avgs4=[avgs4 avg4];
end
average1=[avgs1; 1:11]
average2=[avgs2; 1:11]
average3=[avgs3; 1:11]
average4=[avgs4; 1:11]

tableavg1 = array2table(average1','VariableNames',{'magnitude','position'});
tableavg2 = array2table(average2','VariableNames',{'magnitude','position'});
tableavg3 = array2table(average3','VariableNames',{'magnitude','position'});
tableavg4 = array2table(average4','VariableNames',{'magnitude','position'});
%% 
load('data/data_01-05-2021_13-28-34.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1', 'E5'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope2(filt_datas,Fs);

impactN = 15; % ADJUST number of expected impacts in data

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
avgs1=[];
avgs2=[];
avgs3=[];
avgs4=[];
for i=1:5:15
    avg1= mean(peak_val(1,i:i+4));
    avgs1=[avgs1 avg1];
    avg2= mean(peak_val(2,i:i+4));
    avgs2=[avgs2 avg2];
    avg3= mean(peak_val(3,i:i+4));
    avgs3=[avgs3 avg3];
    avg4= mean(peak_val(4,i:i+4));
    avgs4=[avgs4 avg4];
end

load('data/data_01-05-2021_13-29-59.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1', 'E5'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope2(filt_datas,Fs);

impactN = 15; % ADJUST number of expected impacts in data

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);

for i=1:5:15
    avg1= mean(peak_val(1,i:i+4));
    avgs1=[avgs1 avg1];
    avg2= mean(peak_val(2,i:i+4));
    avgs2=[avgs2 avg2];
    avg3= mean(peak_val(3,i:i+4));
    avgs3=[avgs3 avg3];
    avg4= mean(peak_val(4,i:i+4));
    avgs4=[avgs4 avg4];
end

load('data/data_01-05-2021_13-31-27.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1', 'E5'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope2(filt_datas,Fs);

impactN = 15; % ADJUST number of expected impacts in data

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);

for i=1:5:15
    avg1= mean(peak_val(1,i:i+4));
    avgs1=[avgs1 avg1];
    avg2= mean(peak_val(2,i:i+4));
    avgs2=[avgs2 avg2];
    avg3= mean(peak_val(3,i:i+4));
    avgs3=[avgs3 avg3];
    avg4= mean(peak_val(4,i:i+4));
    avgs4=[avgs4 avg4];
end

load('data/data_01-05-2021_13-32-32.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1', 'E5'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope2(filt_datas,Fs);

impactN = 10; % ADJUST number of expected impacts in data

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);

for i=1:5:10
    avg1= mean(peak_val(1,i:i+4));
    avgs1=[avgs1 avg1];
    avg2= mean(peak_val(2,i:i+4));
    avgs2=[avgs2 avg2];
    avg3= mean(peak_val(3,i:i+4));
    avgs3=[avgs3 avg3];
    avg4= mean(peak_val(4,i:i+4));
    avgs4=[avgs4 avg4];
end
average1_2=[avgs1; 1:11]
average2_2=[avgs2; 1:11]
average3_2=[avgs3; 1:11]
average4_2=[avgs4; 1:11]

