differences1=[];
differences2=[];
differences3=[];

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
steplength = 1;
speeds=[];
for i = 1:length(differences1)
    speed = steplength/differences1(i);
    speeds = [speeds speed];
end
speedno = speeds;

differences1=[];
differences2=[];
differences3=[];
filenames ={'data/data_10-31-2020_17-08-30'};%slight limp
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
filenames ={'data/data_10-31-2020_17-22-35'};%severe limp
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
steplength = 1;
speeds=[];
for i = 1:length(differences1)
    speed = steplength/differences1(i);
    speeds = [speeds speed];
end
speedsevere=speeds;

totalspeeds = [speedno;speedslight;speedsevere];
% err=std(totalspeeds')
% err1=err'
% errs=[];
% for i =1:5
%     errs= [errs err1];
% end
figure
plot(totalspeeds','-o')
title('Speed for each stride')
legend('no limp','slight limp','severe limp')
ylabel('Feet/Sec')
% hold on
% errorbar([1 2 3 4 5; 1 2 3 4 5; 1 2 3 4 5],totalspeeds,errs)
% ylim(1:2)
% totalllength=[0];
% timecom = 0;
% for i = 1:length(differences1)
%     timecom = timecom+differences1(i);
%     timeperstride=[timeperstride timecom];
% end
% stridelengthno=stridelength