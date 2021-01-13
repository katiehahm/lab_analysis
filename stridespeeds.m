function [totalspeeds] = stridespeeds(nolimpfilename,slightlimpfilename, severelimpfilename,steplength)
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

totalspeeds = [speedno;speedslight;speedsevere];

figure
plot(totalspeeds','-o')
title('Speed for each stride')
legend('no limp','slight limp','severe limp')
ylabel('Feet/Sec')