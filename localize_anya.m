%%black soccer ball 31.5" high
data_root_anya = '/Users/anyachase 1/Dropbox (MIT)/Analysis/data/';
load([data_root_anya, 'data_12-10-2020_12-16-21'])
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

% load('data/data_12-10-2020_12-18-11.mat')
% filt_datas = lpf_data(datas);
% loc_names = {'A1', 'A5', 'E1', 'E5'};
% Fs = 12800;
% filt_datas = lpf_data(datas);
% clean_data = clean_envelope2(filt_datas,Fs);
% 
% impactN = 15; % ADJUST number of expected impacts in data
% 
% [onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
% 
% for i=1:5:15
%     avg1= mean(peak_val(1,i:i+4));
%     avgs1=[avgs1 avg1];
%     avg2= mean(peak_val(2,i:i+4));
%     avgs2=[avgs2 avg2];
%     avg3= mean(peak_val(3,i:i+4));
%     avgs3=[avgs3 avg3];
%     avg4= mean(peak_val(4,i:i+4));
%     avgs4=[avgs4 avg4];
% end
% 
% load('data/data_12-10-2020_12-19-54.mat')
% filt_datas = lpf_data(datas);
% loc_names = {'A1', 'A5', 'E1', 'E5'};
% Fs = 12800;
% filt_datas = lpf_data(datas);
% clean_data = clean_envelope2(filt_datas,Fs);
% 
% impactN = 15; % ADJUST number of expected impacts in data
% 
% [onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
% 
% for i=1:5:15
%     avg1= mean(peak_val(1,i:i+4));
%     avgs1=[avgs1 avg1];
%     avg2= mean(peak_val(2,i:i+4));
%     avgs2=[avgs2 avg2];
%     avg3= mean(peak_val(3,i:i+4));
%     avgs3=[avgs3 avg3];
%     avg4= mean(peak_val(4,i:i+4));
%     avgs4=[avgs4 avg4];
% end
% 
% load('data/data_12-10-2020_12-21-17.mat')
% filt_datas = lpf_data(datas);
% loc_names = {'A1', 'A5', 'E1', 'E5'};
% Fs = 12800;
% filt_datas = lpf_data(datas);
% clean_data = clean_envelope2(filt_datas,Fs);
% 
% impactN = 10; % ADJUST number of expected impacts in data
% 
% [onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
% 
% for i=1:5:10
%     avg1= mean(peak_val(1,i:i+4));
%     avgs1=[avgs1 avg1];
%     avg2= mean(peak_val(2,i:i+4));
%     avgs2=[avgs2 avg2];
%     avg3= mean(peak_val(3,i:i+4));
%     avgs3=[avgs3 avg3];
%     avg4= mean(peak_val(4,i:i+4));
%     avgs4=[avgs4 avg4];
% end
% average1=[avgs1; 1:11]
% average2=[avgs2; 1:11]
% average3=[avgs3; 1:11]
% average4=[avgs4; 1:11]
% 
% m=0.98
% PE1=m*9.81*.8
% scalefactor1=PE1./avgs1
% scalefactor2=PE1./avgs2
% scalefactor3=PE1./avgs3
% scalefactor4=PE1./avgs4

% tableavg1 = array2table(average1','VariableNames',{'magnitude','position'});
% tableavg2 = array2table(average2','VariableNames',{'magnitude','position'});
% tableavg3 = array2table(average3','VariableNames',{'magnitude','position'});
% tableavg4 = array2table(average4','VariableNames',{'magnitude','position'});
%% black soccer ball 17" high
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
m=0.98
PE2=m*9.81*0.4318
E2_1 = scalefactor1.*average1_2(1,:)
E2_2 = scalefactor2.*average2_2(1,:)
E2_3 = scalefactor3.*average3_2(1,:)
E2_4 = scalefactor4.*average4_2(1,:)
avgpe=[mean(E2_1) mean(E2_2) mean(E2_3) mean(E2_4)]

%% white soccer ball 31.5" high
load('data/data_01-05-2021_13-36-09.mat')
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

load('data/data_01-05-2021_13-38-22.mat')
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

load('data/data_01-05-2021_13-40-13.mat')
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

load('data/data_01-05-2021_13-41-42.mat')
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
average1_3=[avgs1; 1:11]
average2_3=[avgs2; 1:11]
average3_3=[avgs3; 1:11]
average4_3=[avgs4; 1:11]

tableavg1 = array2table(average1','VariableNames',{'magnitude','position'});
tableavg2 = array2table(average2','VariableNames',{'magnitude','position'});
tableavg3 = array2table(average3','VariableNames',{'magnitude','position'});
tableavg4 = array2table(average4','VariableNames',{'magnitude','position'});

m=0.26
PE3=m*9.81*.8
E3_1 = scalefactor1.*average1_3(1,:)
E3_2 = scalefactor2.*average2_3(1,:)
E3_3 = scalefactor3.*average3_3(1,:)
E3_4 = scalefactor4.*average4_3(1,:)
avgpe3=[mean(E3_1) mean(E3_2) mean(E3_3) mean(E3_4)]
%% white soccer ball trial 2
load('data/data_01-14-2021_19-35-40.mat')
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

load('data/data_01-14-2021_19-39-21.mat')
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

load('data/data_01-14-2021_19-42-34.mat')
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

load('data/data_01-14-2021_19-54-28.mat')
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
average1_3=[avgs1; 1:11]
average2_3=[avgs2; 1:11]
average3_3=[avgs3; 1:11]
average4_3=[avgs4; 1:11]

tableavg1 = array2table(average1','VariableNames',{'magnitude','position'});
tableavg2 = array2table(average2','VariableNames',{'magnitude','position'});
tableavg3 = array2table(average3','VariableNames',{'magnitude','position'});
tableavg4 = array2table(average4','VariableNames',{'magnitude','position'});

m=0.26
PE3=m*9.81*0.635
E3_1 = scalefactor1.*average1_3(1,:)
E3_2 = scalefactor2.*average2_3(1,:)
E3_3 = scalefactor3.*average3_3(1,:)
E3_4 = scalefactor4.*average4_3(1,:)
avgpe3=[mean(E3_1) mean(E3_2) mean(E3_3) mean(E3_4)]

%% 
m=1
PE1=m*9.81*.8
% sensor1 = [average1(1,:);average1_2(1,:);average1_3]
% sensor2 = [average2(1,:);average2_2(1,:);average2_3]
% sensor3 = [average3(1,:);average3_2(1,:);average3_3]
% sensor4 = [average4(1,:);average4_2(1,:);average4_3]
% tableavg1 = array2table(sensor1,'VariableNames',{'1','2','3','4','5','6','7','8','9','10','11'});
% tableavg2 = array2table(sensor2,'VariableNames',{'1','2','3','4','5','6','7','8','9','10','11'});
% tableavg3 = array2table(sensor3,'VariableNames',{'1','2','3','4','5','6','7','8','9','10','11'});
% tableavg4 = array2table(sensor4,'VariableNames',{'1','2','3','4','5','6','7','8','9','10','11'});

%% 
[X,Y]=meshgrid(1:7);
figure; hold on;
plot(X,Y,'k');
plot(Y,X,'k');axis off
I=(rand(7));
surface(I);
% h=linspace(0.5,1,64);
% h=[h',h',h'];
% set(gcf,'Colormap',h);
% N=36:-1:1;q=1;
% x=linspace(1.5,10.5,10);
% y=linspace(1.5,10.5,10);
% for n=1:6
%     for p=1:6
%         text(y(n)-.5,x(p),num2str(N(q)),'FontWeight','bold');
%         q=q+1;
%     end
% end
%% %% grid mapping function
filenames={'data_01-05-2021_13-28-34.mat', 'data_01-05-2021_13-29-59.mat', 'data_01-05-2021_13-31-27.mat', 'data_01-05-2021_13-32-32.mat'}
impacts=[15 15 15 10]
bounceN=5
height=0.4318
mass=.98
peaks=[]
[error1, error2, error3, error4] = grid_function(filenames,impacts,bounceN,mass,height)
