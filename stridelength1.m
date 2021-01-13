%length of stride
differences1=[];
differences2=[];
differences3=[];
% filenames = {'data/data_10-19-2020_15-49.mat', 'data/data_10-19-2020_15-56.mat','data/data_10-19-2020_15-57.mat','data/data_10-19-2020_15-58.mat','data/data_10-19-2020_15-59.mat','data/data_10-19-2020_16-00.mat',...
%     'data/data_10-19-2020_16-01.mat','data/data_10-19-2020_16-02.mat','data/data_10-19-2020_16-03.mat','data/data_10-19-2020_16-04.mat','data/data_10-19-2020_16-05.mat',...
%     'data/data_10-19-2020_16-06.mat','data/data_10-19-2020_16-07.mat','data/data_10-19-2020_16-08.mat','data/data_10-19-2020_16-09.mat','data/data_10-19-2020_16-10.mat','data/data_10-19-2020_16-11.mat',...
%     'data/data_10-19-2020_16-12.mat','data/data_10-19-2020_16-13.mat','data/data_10-19-2020_16-14.mat','data/data_10-19-2020_16-15.mat','data/data_10-19-2020_16-16.mat','data/data_10-19-2020_16-17.mat','data/data_10-19-2020_16-18.mat',...
%     'data/data_10-19-2020_16-19.mat','data/data_10-19-2020_16-20.mat','data/data_10-19-2020_16-21.mat','data/data_10-19-2020_16-22.mat','data/data_10-19-2020_16-23.mat','data/data_10-19-2020_16-24.mat',...
%     'data/data_10-19-2020_16-25.mat','data/data_10-19-2020_16-26.mat','data/data_10-19-2020_16-27.mat','data/data_10-19-2020_16-28.mat','data/data_10-19-2020_16-29.mat','data/data_10-19-2020_16-30.mat',...
%     'data/data_10-19-2020_16-31.mat','data/data_10-19-2020_16-32.mat','data/data_10-19-2020_16-33.mat','data/data_10-19-2020_16-34.mat','data/data_10-19-2020_16-35.mat','data/data_10-19-2020_16-36.mat',...
%     'data/data_10-19-2020_16-37.mat','data/data_10-19-2020_16-38.mat','data/data_10-19-2020_16-39.mat','data/data_10-19-2020_16-40.mat','data/data_10-19-2020_16-41.mat','data/data_10-19-2020_16-42.mat','data/data_10-19-2020_16-43.mat','data/data_10-19-2020_16-44.mat'}; 
filenames ={'data/data_10-19-2020_15-49.mat'};%no limp
for kk = 1:numel(filenames)
    load(filenames{kk});
    filt_datas = lpf_data(datas);
    loc_names = {'A1', 'A5', 'E1'};
    Fs = 12800;
    clean_data = clean_envelope(filt_datas,Fs);
    impactN = 6; % variable input value

    [onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
    [diffsinlength1,diffsinlength2,diffsinlength3] = cadence(peak_idx,impactN,Fs);
    differences1= [differences1;diffsinlength1];
    differences2= [differences2;diffsinlength2];
    differences3= [differences3;diffsinlength3];
end
totallength = 6;
totaltime=sum(differences1);
speed = totallength/totaltime;
location= 0;
stridelength=[0];
for i = 1:length(differences1)
    location = location+differences1(i)*speed;
    stridelength = [stridelength location];
end
timeperstride=[0];
timecom = 0;
for i = 1:length(differences1)
    timecom = timecom+differences1(i);
    timeperstride=[timeperstride timecom];
end
stridelengthno=stridelength
% stridelength
% timeperstride
% stridelength./timeperstride
differences1=[];
differences2=[];
differences3=[];
figure
filenames ={'data/data_10-31-2020_17-22-35'}; %severe limp
for kk = 1:numel(filenames)
    load(filenames{kk});
    filt_datas = lpf_data(datas);
    loc_names = {'A1', 'A5', 'E1'};
    Fs = 12800;
    clean_data = clean_envelope(filt_datas,Fs);
    impactN = 6; % variable input value

    [onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
    [diffsinlength1,diffsinlength2,diffsinlength3] = cadence(peak_idx,impactN,Fs);
    differences1= [differences1;diffsinlength1];
    differences2= [differences2;diffsinlength2];
    differences3= [differences3;diffsinlength3];
end
totallength = 6;
totaltime=sum(differences1);
speed = totallength/totaltime;
location= 0;
stridelength=[0];
for i = 1:length(differences1)
    location = location+differences1(i)*speed;
    stridelength = [stridelength location];
end
timeperstride=[0];
timecom = 0;
for i = 1:length(differences1)
    timecom = timecom+differences1(i);
    timeperstride=[timeperstride timecom];
end
stridelengthsevere=stridelength
differences1=[];
differences2=[];
differences3=[];
figure
filenames ={'data/data_10-31-2020_17-08-30'}; %slight limp
for kk = 1:numel(filenames)
    load(filenames{kk});
    filt_datas = lpf_data(datas);
    loc_names = {'A1', 'A5', 'E1'};
    Fs = 12800;
    clean_data = clean_envelope(filt_datas,Fs);
    impactN = 6; % variable input value

    [onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
    [diffsinlength1,diffsinlength2,diffsinlength3] = cadence(peak_idx,impactN,Fs);
    differences1= [differences1;diffsinlength1];
    differences2= [differences2;diffsinlength2];
    differences3= [differences3;diffsinlength3];
end
totallength = 6;
totaltime=sum(differences1);
speed = totallength/totaltime;
location= 0;
stridelength=[0];
for i = 1:length(differences1)
    location = location+differences1(i)*speed;
    stridelength = [stridelength location];
end
timeperstride=[0]
timecom = 0;
for i = 1:length(differences1)
    timecom = timecom+differences1(i);
    timeperstride=[timeperstride timecom];
end
stridelengthslight=stridelength

strides=[stridelengthno; stridelengthslight; stridelengthsevere];
figure
bar(strides')
legend('nolimp','slight limp','severelimp')
title('Distance Travelled')
ylabel('length in feet')