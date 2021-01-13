function [diff_foot1, diff_foot2] = hist_stridesamefoot(filenames,impactN)
diff_foot1=[];
diff_foot2=[];
for kk = 1:numel(filenames)
    foot1=[];
    foot2=[];
    load(filenames{kk});
    filt_datas = lpf_data(datas);
    loc_names = {'A1', 'A5', 'E1'};
    Fs = 12800;
    clean_data = clean_envelope(filt_datas,Fs);

    [onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names,false);
    peak_idx
    sortpeakidx= [sort(peak_idx(1,:));sort(peak_idx(2,:));sort(peak_idx(3,:))]
    for i=1:2:impactN
        foot1=[foot1 sortpeakidx(:,i)];
    end
    for i=2:2:impactN
        foot2=[foot2 sortpeakidx(:,i)];
    end
    for i=2:length(foot1(1,:))
        diff_foot1=[diff_foot1 foot1(:,i)-foot1(:,i-1)];
    end
    for i=2:length(foot2(2,:))
        diff_foot2=[diff_foot2 foot2(:,i)-foot2(:,i-1)];
    end    
end
diff_foot1=diff_foot1./Fs;
diff_foot2=diff_foot2./Fs;
bin = round(1+3.22*log(numel(diff_foot1(1,:))));
figure
hf=histfit(diff_foot1(1,:),bin,'kernel');
figure
histfit(diff_foot1(1,:),bin)
hold on
x=get(hf(2),'XData'); 
y=get(hf(2),'YData');
plot(x,y,'Color','g','LineWidth', 2)
title('Stride Time Sensor 1 Foot1')
xlabel('Stride Time')
ylabel('Occurances')
mu1=mean(diff_foot1(1,:));
sigma1=std(diff_foot1(1,:));
hold on
line([mu1, mu1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu1 + sigma1, mu1 + sigma1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu1 - sigma1, mu1 - sigma1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
hold off

 bin = round(1+3.22*log(numel(diff_foot2(1,:))));
figure
hf=histfit(diff_foot2(1,:),bin,'kernel');
figure
histfit(diff_foot2(1,:),bin)
hold on
x=get(hf(2),'XData'); 
y=get(hf(2),'YData');
plot(x,y,'Color','g','LineWidth', 2)
title('Stride Time Sensor 1 Foot2')
xlabel('Stride Time')
ylabel('Occurances')
mu1=mean(diff_foot2(1,:));
sigma1=std(diff_foot2(1,:));
hold on
line([mu1, mu1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu1 + sigma1, mu1 + sigma1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu1 - sigma1, mu1 - sigma1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
hold off

bin = round(1+3.22*log(numel(diff_foot1(2,:))));
figure
hf=histfit(diff_foot1(2,:),bin,'kernel');
figure
histfit(diff_foot1(2,:),bin)
hold on
x=get(hf(2),'XData'); 
y=get(hf(2),'YData');
plot(x,y,'Color','g','LineWidth', 2)
title('Stride Time Sensor 2 Foot1')
xlabel('Stride Time')
ylabel('Occurances')
mu1=mean(diff_foot1(2,:));
sigma1=std(diff_foot1(2,:));
hold on
line([mu1, mu1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu1 + sigma1, mu1 + sigma1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu1 - sigma1, mu1 - sigma1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
hold off

bin = round(1+3.22*log(numel(diff_foot2(2,:))));
figure
hf=histfit(diff_foot2(2,:),bin,'kernel');
figure
histfit(diff_foot2(2,:),bin)
hold on
x=get(hf(2),'XData'); 
y=get(hf(2),'YData');
plot(x,y,'Color','g','LineWidth', 2)
title('Stride Time Sensor 2 Foot2')
xlabel('Stride Time')
ylabel('Occurances')
mu1=mean(diff_foot2(2,:));
sigma1=std(diff_foot2(2,:));
hold on
line([mu1, mu1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu1 + sigma1, mu1 + sigma1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu1 - sigma1, mu1 - sigma1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
hold off

bin = round(1+3.22*log(numel(diff_foot1(3,:))));
figure
hf=histfit(diff_foot1(3,:),bin,'kernel');
figure
histfit(diff_foot1(3,:),bin)
hold on
x=get(hf(2),'XData'); 
y=get(hf(2),'YData');
plot(x,y,'Color','g','LineWidth', 2)
title('Stride Time Sensor 3 Foot1')
xlabel('Stride Time')
ylabel('Occurances')
mu1=mean(diff_foot1(3,:));
sigma1=std(diff_foot1(3,:));
hold on
line([mu1, mu1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu1 + sigma1, mu1 + sigma1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu1 - sigma1, mu1 - sigma1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
hold off

bin = round(1+3.22*log(numel(diff_foot2(3,:))));
figure
hf=histfit(diff_foot2(3,:),bin,'kernel');
figure
histfit(diff_foot2(3,:),bin)
hold on
x=get(hf(2),'XData'); 
y=get(hf(2),'YData');
plot(x,y,'Color','g','LineWidth', 2)
title('Stride Time Sensor 3 Foot2')
xlabel('Stride Time')
ylabel('Occurances')
mu1=mean(diff_foot2(3,:));
sigma1=std(diff_foot2(3,:));
hold on
line([mu1, mu1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu1 + sigma1, mu1 + sigma1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
line([mu1 - sigma1, mu1 - sigma1], ylim, 'Color', 'c', 'LineWidth', 0.5); 
hold off
end