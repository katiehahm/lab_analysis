%%
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj';
subj = '1'; % number of subject
intervention = 'normal2'; % normal1 slow insole weight count normal2
filename = [data_root_katie, subj, '_', intervention];

load(filename)
Fs = 12800;
differences = zeros(1,length(arrival_idx)-1);
for i = 1:(length(arrival_idx)-1)
    curr = min(arrival_idx(i+1,:));
    prev = min(arrival_idx(i,:));
    differences(i) = (curr - prev)/Fs;
end

bin = round(1+3.22*log(numel(differences)));
% figure
% hf=histfit(differences,bin,'kernel');
figure
histfit(differences,bin)
hold on
x=get(hf(2),'XData'); 
y=get(hf(2),'YData');
plot(x,y,'Color','k','LineWidth', 2)
titleroot = 'Step time distribution ';
title([titleroot, intervention])
xlabel('Step Time')
ylabel('Occurances')
mu=mean(differences);
sigma=std(differences);
hold on
line([mu, mu], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu + sigma, mu + sigma], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu - sigma, mu - sigma], ylim, 'Color', 'c', 'LineWidth', 0.5); 
legend('Histogram','Normal','Fitted','1 std')

mu
sigma
[BF, BC] = bimodalitycoeff(differences)

%% old stuff from other folder
differences1=[];
differences2=[];
differences3=[];
for kk = 1:numel(filenames)
    load(filenames{kk});
    filt_datas = lpf_data(datas);
    loc_names = {'A1', 'A5', 'E1'};
    Fs = 12800;
    clean_data = clean_envelope(filt_datas,Fs);

    [onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names,false);
    [diffsinlength1,diffsinlength2,diffsinlength3] = cadence(peak_idx,impactN,Fs);
    differences1= [differences1;diffsinlength1];
    differences2= [differences2;diffsinlength2];
    differences3= [differences3;diffsinlength3];
end

bin = round(1+3.22*log(numel(differences1)));
figure
hf=histfit(differences1,bin,'kernel');
figure
histfit(differences1,bin)
hold on
x=get(hf(2),'XData'); 
y=get(hf(2),'YData');
plot(x,y,'Color','g','LineWidth', 2)
title('Stride Time Sensor 1')
xlabel('Stride Time')
ylabel('Occurances')
mu=mean(differences1);
sigma=std(differences1);
hold on
line([mu, mu], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu + sigma, mu + sigma], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu - sigma, mu - sigma], ylim, 'Color', 'c', 'LineWidth', 0.5); 
hold off

bin = round(1+3.22*log(numel(differences2)));
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
mu2=mean(differences2);
sigma2=std(differences2);
hold on
line([mu2, mu2], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu2 + sigma2, mu2 + sigma2], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu2 - sigma2, mu2 - sigma2], ylim, 'Color', 'c', 'LineWidth', 0.5); 
bin = round(1+3.22*log(numel(differences3)));
hold off

bin = round(1+3.22*log(numel(differences3)));
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
mu3=mean(differences3);
sigma3=std(differences3);
hold on
line([mu3, mu3], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu3 + sigma3, mu3 + sigma3], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu3 - sigma3, mu3 - sigma3], ylim, 'Color', 'c', 'LineWidth', 0.5); 
