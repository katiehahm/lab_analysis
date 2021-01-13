nolimpfilename ={'data/data_10-19-2020_15-49.mat'};%no limp
slightlimpfilename ={'data/data_10-31-2020_17-08-30'}; %slight limp
severelimpfilename ={'data/data_10-31-2020_17-22-35'}; %severe limp
steplength = 1;%consistent step length

differences1=[];
differences2=[];
differences3=[];
for kk = 1:numel(nolimpfilename)
    load(nolimpfilename{kk});
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
steplength=1;
speeds=[];
for i = 1:length(differences1)
    speed = steplength/differences1(i);
    speeds = [speeds speed];
end
speedno = speeds;

differences1=[];
differences2=[];
differences3=[];
for kk = 1:numel(slightlimpfilename)
    load(slightlimpfilename{kk});
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
steplength = 1;
speeds=[];
for i = 1:length(differences1)
    speed = steplength/differences1(i);
    speeds = [speeds speed];
end
speedslight=speeds;

differences1=[];
differences2=[];
differences3=[];
for kk = 1:numel(severelimpfilename)
    load(severelimpfilename{kk});
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
steplength = 1;
speeds=[];
for i = 1:length(differences1)
    speed = steplength/differences1(i);
    speeds = [speeds speed];
end
speedsevere=speeds;

speedconstant= {'data/data_11-06-2020_11-58-08'};
differences1=[];
differences2=[];
differences3=[];
for kk = 1:numel(speedconstant)
    load(speedconstant{kk});
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

steplength = [18/12 18/12 18/12 18/12 18/12]; %[14/12 22/12 14/12 22/12 14/12];
speedconst=differences1./steplength';

steplength=[14/12 22/12 14/12 22/12 14/12];
limpspeed = steplength'./differences1;
% %totalspeeds = [speedno;speedslight;speedsevere;speedconst'];
% totalspeeds = [speedconst'; limpspeed'];
% figure
% plot(totalspeeds','-o')
% title('Speed for each stride')
% %legend('no limp','slight limp','severe limp','constantspeed')
% legend('idealwalkspeed', 'actualwalkspeed')
% ylabel('Feet/Sec')

figure
differences1=[];
differences2=[];
differences3=[];
for kk = 1:numel(nolimpfilename)
    load(nolimpfilename{kk});
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
totallength=18/12*5
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

differences1=[];
differences2=[];
differences3=[];
figure

for kk = 1:numel(severelimpfilename)
    load(severelimpfilename{kk});
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
totallength = 18/12*5;
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
for kk = 1:numel(slightlimpfilename)
    load(slightlimpfilename{kk});
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
totallength = 18/12*5;
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
stridelengthslight=stridelength

differences1=[];
differences2=[];
differences3=[];
for kk = 1:numel(speedconstant)
    load(speedconstant{kk});
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

location= 0;
stridelength=[0];
for i = 1:length(differences1)
    location = location+differences1(i)*limpspeed(i);
    stridelength = [stridelength location];
end
newlimp=stridelength

figure
strides=[stridelengthno; stridelengthslight; stridelengthsevere;newlimp];
figure
bar(strides')
legend('nolimp','slight limp','severelimp','newlimp')
title('Distance Travelled')
ylabel('length in feet')
