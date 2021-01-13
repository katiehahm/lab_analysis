function [differences1, differences2, differences3, mu1, mu2,mu3,sigma1,sigma2,sigma3] = hist_stridetime(filenames,impactN)
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
mu1=mean(differences1);
sigma1=std(differences1);
hold on
line([mu1, mu1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu1 + sigma1, mu1 + sigma1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu1 - sigma1, mu1 - sigma1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
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
end