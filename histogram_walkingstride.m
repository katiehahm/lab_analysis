filenames = {'data/data_10-19-2020_15-49.mat', 'data/data_10-19-2020_15-56.mat','data/data_10-19-2020_15-57.mat','data/data_10-19-2020_15-58.mat','data/data_10-19-2020_15-59.mat','data/data_10-19-2020_16-00.mat',...
    'data/data_10-19-2020_16-01.mat','data/data_10-19-2020_16-02.mat','data/data_10-19-2020_16-03.mat','data/data_10-19-2020_16-04.mat','data/data_10-19-2020_16-05.mat',...
    'data/data_10-19-2020_16-06.mat','data/data_10-19-2020_16-07.mat','data/data_10-19-2020_16-08.mat','data/data_10-19-2020_16-09.mat','data/data_10-19-2020_16-10.mat','data/data_10-19-2020_16-11.mat',...
    'data/data_10-19-2020_16-12.mat','data/data_10-19-2020_16-13.mat','data/data_10-19-2020_16-14.mat','data/data_10-19-2020_16-15.mat','data/data_10-19-2020_16-16.mat','data/data_10-19-2020_16-17.mat','data/data_10-19-2020_16-18.mat',...
    'data/data_10-19-2020_16-19.mat','data/data_10-19-2020_16-20.mat','data/data_10-19-2020_16-21.mat','data/data_10-19-2020_16-22.mat','data/data_10-19-2020_16-23.mat','data/data_10-19-2020_16-24.mat',...
    'data/data_10-19-2020_16-25.mat','data/data_10-19-2020_16-26.mat','data/data_10-19-2020_16-27.mat','data/data_10-19-2020_16-28.mat','data/data_10-19-2020_16-29.mat','data/data_10-19-2020_16-30.mat',...
    'data/data_10-19-2020_16-31.mat','data/data_10-19-2020_16-32.mat','data/data_10-19-2020_16-33.mat','data/data_10-19-2020_16-34.mat','data/data_10-19-2020_16-35.mat','data/data_10-19-2020_16-36.mat',...
    'data/data_10-19-2020_16-37.mat','data/data_10-19-2020_16-38.mat','data/data_10-19-2020_16-39.mat','data/data_10-19-2020_16-40.mat','data/data_10-19-2020_16-41.mat','data/data_10-19-2020_16-42.mat','data/data_10-19-2020_16-43.mat','data/data_10-19-2020_16-44.mat'}; 
% filenames = {'data/data_10-19-2020_16-50.mat','data/data_10-19-2020_16-51.mat'};
differences1=[];
differences2=[];
differences3=[];
%load('data/data_10-19-2020_15-49.mat') %data_10-19-2020_15-49 -- 16-44 walk BQ,CQ,DQ,EQ,ER,ES in that order. 6 steps.
%   load('data/data_10-19-2020_15-49.mat');  
%   filt_datas = lpf_data(datas);
%     loc_names = {'A1', 'A5', 'E1'};
%     Fs = 12800;
%     clean_data = clean_envelope(filt_datas,Fs);
%     impactN = 6; % variable input value
% 
%     [onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
%     [cadencetotal,avgstridelengths,stddevstridelength,diffsinlength1,diffsinlength2,diffsinlength3] = cadence(peak_idx,impactN,Fs);
%     differences1= diffsinlength1;
%     differences2= diffsinlength2;
%     differences3= diffsinlength3;
% 
% datastr = 'data/data_10-19-2020_15-';
% strmat = '.mat';
% 
% for i = 56:59
%     loadstr= strcat(strcat(datastr, int2str(i)) , strmat);
%     load(loadstr)
%      filt_datas = lpf_data(datas);
%     loc_names = {'A1', 'A5', 'E1'};
%     Fs = 12800;
%     clean_data = clean_envelope(filt_datas,Fs);
%     impactN = 6; % variable input value
% 
%     [onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
%     [cadencetotal,avgstridelengths,stddevstridelength,diffsinlength1,diffsinlength2,diffsinlength3] = cadence(peak_idx,impactN,Fs);
%     differences1= [differences1;diffsinlength1];
%     differences2= [differences2;diffsinlength2];
%     differences3= [differences3;diffsinlength3];
% end
% % 
% datastr = 'data/data_10-19-2020_16-0';
% strmat = '.mat';
% 
% for i = 0:9
%     loadstr= strcat(strcat(datastr, int2str(i)) , strmat);
%     load(loadstr)
%      filt_datas = lpf_data(datas);
%     loc_names = {'A1', 'A5', 'E1'};
%     Fs = 12800;
%     clean_data = clean_envelope(filt_datas,Fs);
%     impactN = 6; % variable input value
% 
%     [onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
%     [cadencetotal,avgstridelengths,stddevstridelength,diffsinlength1,diffsinlength2,diffsinlength3] = cadence(peak_idx,impactN,Fs);
%     differences1= [differences1;diffsinlength1];
%     differences2= [differences2;diffsinlength2];
%     differences3= [differences3;diffsinlength3];
% end
% datastr = 'data/data_10-19-2020_16-';
% strmat = '.mat';
% 
% for i = 10:44
%     loadstr= strcat(strcat(datastr, int2str(i)) , strmat);
%     load(loadstr)
%      filt_datas = lpf_data(datas);
%     loc_names = {'A1', 'A5', 'E1'};
%     Fs = 12800;
%     clean_data = clean_envelope(filt_datas,Fs);
%     impactN = 6; % variable input value
% 
%     [onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
%     [cadencetotal,avgstridelengths,stddevstridelength,diffsinlength1,diffsinlength2,diffsinlength3] = cadence(peak_idx,impactN,Fs);
%     differences1= [differences1;diffsinlength1];
%     differences2= [differences2;diffsinlength2];
%     differences3= [differences3;diffsinlength3];
% end
% 
% %slight limp
%  filenames = {'data/data_10-31-2020_17-08-30','data/data_10-31-2020_17-09-11','data/data_10-31-2020_17-10-07','data/data_10-31-2020_17-10-43','data/data_10-31-2020_17-11-18','data/data_10-31-2020_17-12-30',...
%     'data/data_10-31-2020_17-13-06','data/data_10-31-2020_17-13-45','data/data_10-31-2020_17-14-14','data/data_10-31-2020_17-14-53','data/data_10-31-2020_17-15-19','data/data_10-31-2020_17-15-55',...
%     'data/data_10-31-2020_17-16-24','data/data_10-31-2020_17-17-12','data/data_10-31-2020_17-17-43','data/data_10-31-2020_17-18-11','data/data_10-31-2020_17-18-48','data/data_10-31-2020_17-19-11',...
%     'data/data_10-31-2020_17-19-36','data/data_10-31-2020_17-20-04'};

%severe limp
% filenames ={'data/data_10-31-2020_17-22-35','data/data_10-31-2020_17-23-01','data/data_10-31-2020_17-23-26','data/data_10-31-2020_17-23-48','data/data_10-31-2020_17-24-16','data/data_10-31-2020_17-24-38',...
%     'data/data_10-31-2020_17-25-04','data/data_10-31-2020_17-25-29','data/data_10-31-2020_17-25-53','data/data_10-31-2020_17-26-23','data/data_10-31-2020_17-26-48','data/data_10-31-2020_17-27-16',...
%     'data/data_10-31-2020_17-27-45','data/data_10-31-2020_17-28-13','data/data_10-31-2020_17-28-45','data/data_10-31-2020_17-29-19','data/data_10-31-2020_17-29-42','data/data_10-31-2020_17-30-08',...
%     'data/data_10-31-2020_17-30-36','data/data_10-31-2020_17-31-10'};
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
% % figure
% % histogram(strides(:,1),6)
% % title('Average Stride Length Sensor 1')
% % xlabel('Average Stride Length')
% % ylabel('Occurances')
% 
% 

%
% mu=mean(differences1)
% sigma=std(differences1);
% idx= max(find(diffsort<=mu))
% diffsort = sort(differences1);
% bin = round(1+3.22*log(numel(differences1)));
% figure
% hf=histfit(differences1,bin,'kernel');
% figure
% histfit(diffsort(1:idx),bin/2);
% hold on
% histfit(diffsort(idx+1:end),bin/2);
% x=get(hf(2),'XData'); 
% y=get(hf(2),'YData');
% plot(x,y,'Color','g','LineWidth', 2)
% title('Stride Time Sensor 1')
% xlabel('Stride Time')
% ylabel('Occurances')
% hold on
% line([mu, mu], ylim, 'Color', 'c', 'LineWidth', 0.5); 
% line([mu + sigma, mu + sigma], ylim, 'Color', 'c', 'LineWidth', 0.5); 
% line([mu - sigma, mu - sigma], ylim, 'Color', 'c', 'LineWidth', 0.5); 
% 
% diffsort = sort(differences2);
% bin = round(1+3.22*log(numel(differences2)));
% figure
% hf=histfit(differences2,bin,'kernel');
% figure
% histfit(diffsort(1:end/2),bin/2);
% hold on
% histfit(diffsort(end/2:end),bin/2);
% x=get(hf(2),'XData'); 
% y=get(hf(2),'YData');
% plot(x,y,'Color','g','LineWidth', 2)
% title('Stride Time Sensor 2')
% xlabel('Stride Time')
% ylabel('Occurances')
% mu=mean(differences2);
% sigma=std(differences2);
% hold on
% line([mu, mu], ylim, 'Color', 'c', 'LineWidth', 0.5); 
% line([mu + sigma, mu + sigma], ylim, 'Color', 'c', 'LineWidth', 0.5); 
% line([mu - sigma, mu - sigma], ylim, 'Color', 'c', 'LineWidth', 0.5);
% 
% diffsort = sort(differences3);
% bin = round(1+3.22*log(numel(differences3)));
% figure
% hf=histfit(differences3,bin,'kernel');
% figure
% histfit(diffsort(1:end/2),bin/2);
% hold on
% histfit(diffsort(end/2:end),bin/2);
% x=get(hf(2),'XData'); 
% y=get(hf(2),'YData');
% plot(x,y,'Color','g','LineWidth', 2)
% title('Stride Time Sensor 3')
% xlabel('Stride Time')
% ylabel('Occurances')
% mu=mean(differences3);
% sigma=std(differences1);
% hold on
% line([mu, mu], ylim, 'Color', 'c', 'LineWidth', 0.5); 
% line([mu + sigma, mu + sigma], ylim, 'Color', 'c', 'LineWidth', 0.5); 
% line([mu - sigma, mu - sigma], ylim, 'Color', 'c', 'LineWidth', 0.5);


bin = round(1+3.22*log(numel(differences2)));
hold off
figure
hf=histfit(differences2,bin,'kernel');
figure
histfit(differences2,bin)
hold on
x=get(hf(2),'XData'); 
y=get(hf(2),'YData');
plot(x,y,'Color','g','LineWidth', 2)
title('Stride Time Sensor 2')
xlabel('Stride Time')
ylabel('Occurances')
mu=mean(differences2);
sigma=std(differences2);
hold on
line([mu, mu], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu + sigma, mu + sigma], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu - sigma, mu - sigma], ylim, 'Color', 'c', 'LineWidth', 0.5); 
bin = round(1+3.22*log(numel(differences3)));
hold off

figure
hf=histfit(differences3,bin,'kernel');
figure
histfit(differences3,bin)
hold on
x=get(hf(2),'XData'); 
y=get(hf(2),'YData');
plot(x,y,'Color','g','LineWidth', 2)
title('Stride Time Sensor 3')
xlabel('Stride Time')
ylabel('Occurances')
mu=mean(differences3);
sigma=std(differences3);
hold on
line([mu, mu], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu + sigma, mu + sigma], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu - sigma, mu - sigma], ylim, 'Color', 'c', 'LineWidth', 0.5); 