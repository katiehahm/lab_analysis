function [peaks1, peaks2, peaks3,mu1,sigma1, mu2, sigma2,mu3,sigma3] = hist_stridemag(filenames,impactN)
peaks1=[];
peaks2=[];
peaks3=[];
for kk = 1:numel(filenames)
    load(filenames{kk});
    filt_datas = lpf_data(datas);
    loc_names = {'A1', 'A5', 'E1'};
    Fs = 12800;
    clean_data = clean_envelope(filt_datas,Fs);

    [onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names,false);
   peaks1 = [peaks1 peak_val(1,:)];
   peaks2 = [peaks2 peak_val(2,:)];
   peaks3 = [peaks3 peak_val(3,:)];
end

bin = round(1+3.22*log(numel(peaks1)));
figure
hf=histfit(peaks1,bin,'kernel');
figure
histfit(peaks1,bin)
hold on
x=get(hf(2),'XData'); 
y=get(hf(2),'YData');
plot(x,y,'Color','g','LineWidth', 2)
title('Stride Magnitude Sensor 1')
xlabel('Stride Magnitude Sensor 1 (V)')
ylabel('Occurances')
mu1=mean(peaks1);
sigma1=std(peaks1);
hold on
line([mu1, mu1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu1 + sigma1, mu1 + sigma1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu1 - sigma1, mu1 - sigma1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
hold off

bin = round(1+3.22*log(numel(peaks2)));
figure
hf=histfit(peaks2,bin,'kernel');
figure
histfit(peaks2,bin)
hold on
x=get(hf(2),'XData'); 
y=get(hf(2),'YData');
plot(x,y,'Color','g','LineWidth', 2)
title('Stride Magnitude Sensor 2')
xlabel('Stride Magnitude Sensor 2 (V)')
ylabel('Occurances')
mu2=mean(peaks2);
sigma2=std(peaks2);
hold on
line([mu2, mu2], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu2 + sigma2, mu2 + sigma2], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu2 - sigma2, mu2 - sigma2], ylim, 'Color', 'c', 'LineWidth', 0.5); 
hold off

bin = round(1+3.22*log(numel(peaks3)));
figure
hf=histfit(peaks3,bin,'kernel');
figure
histfit(peaks3,bin)
hold on
x=get(hf(2),'XData'); 
y=get(hf(2),'YData');
plot(x,y,'Color','g','LineWidth', 2)
title('Stride Magnitude Sensor 3')
xlabel('Stride Magnitude Sensor 3 (V)')
ylabel('Occurances')
mu3=mean(peaks3);
sigma3=std(peaks3);
hold on
line([mu3, mu3], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu3 + sigma3, mu3 + sigma3], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu3 - sigma3, mu3 - sigma3], ylim, 'Color', 'c', 'LineWidth', 0.5); 
end