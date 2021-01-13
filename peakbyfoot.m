function [avgpeaks1,avgpeaks2,avgpeaktable]=peakbyfoot(filenames,impactN)
avgpeaks1=[];
avgpeaks2=[];
avgdist1 = [];
avgdist2 = [];
X=[];
for kk = 1:numel(filenames)
    load(filenames{kk});
    filt_datas = lpf_data(datas);
    loc_names = {'A1', 'A5', 'E1'};
    Fs = 12800;
    clean_data = clean_envelope(filt_datas,Fs);
  

    [onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
    [diffsinlength1,diffsinlength2,diffsinlength3] = cadence(peak_idx,impactN,Fs);
    differences1= diffsinlength1'
    differences2= diffsinlength2'
    differences3= diffsinlength3'

    for i = 1:2:impactN-1
       avgdist1 = [avgdist1 mean([differences1(i),differences2(i),differences3(i)])];
   end
  for i=2:2:impactN-1
      avgdist2 = [avgdist2 mean([differences1(i),differences2(i),differences3(i)])];
  end

   peaks1 = peak_val(1,:);
   peaks2 = peak_val(2,:);
   peaks3 = peak_val(3,:);
   for i = 1:2:impactN
       avgs = mean([peaks1(i),peaks2(i),peaks3(i)]);
       avgpeaks1= [avgpeaks1 avgs];
   end
  for i=2:2:impactN-1
      avgpeaks2= [avgpeaks2 mean([peaks1(i),peaks2(i),peaks3(i)])];
  end
end

% peakfoot=[avgpeaks1' zeros(numel(avgpeak1))' avgpeaks2' ones(numel(avgpeak2))']
% avgpeaktable = array2table(peakfoot,'VariableNames',{'Sensor_1','Sensor_2','Sensor3','foot'});
% % figure
% % scatter(avgpeaks1,avgpeaks2)
% % title('Peak trends')
% % xlabel('Foot1 Peak Mag')
% % ylabel('Foot2 Peak Mag')
% % ylim([0 0.35])
% % xlim([0 0.35])
% % hold off
% % figure
% % scatter(avgpeaks1,avgpeaks1)
% % hold on
% % scatter(avgpeaks2,avgpeaks2)
% % legend('foot1','foot2')
% % title('Peak by foot')
% % xlabel('Avg Peak Mag')
% % ylabel('Avg Peak Mag')
% % ylim([0 0.35])
% % xlim([0 0.35])
dist=[avgdist1'; avgdist2'];
peak=[avgpeaks1'; avgpeaks2'];
foot1=ones(numel(avgdist1),1);
foot2=2.*ones(numel(avgdist2),1);
foot=[foot1;foot2];
 peakdistfoot = [peak dist foot];
 avgpeaktable = array2table(peakdistfoot,'VariableNames',{'Avg Peaks','Avg Dist','foot'});
X=[avgpeaks1 avgpeaks2;avgpeaks1 avgpeaks2] %%Gausian Distribution below
 X = [peak dist]
figure
scatter(X(:,1),X(:,2));
% scatter(avgdist1,avgpeaks1)
% hold on
% scatter(avgdist2,avgpeaks2)
title('Stride Dist vs Peak Mag')
% legend('foot1','foot2')
xlabel('Avg Stride Dist')
ylabel('Avg Peak Mag')

GMModel = fitgmdist(X,2);
figure
y = [zeros(numel(avgdist1),1);ones(numel(avgdist2),1)]
h = gscatter(X(:,1),X(:,2),y);
hold on
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
g = gca;
fcontour(gmPDF,[g.XLim g.YLim])
title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
legend(h,'Foot 1','Foot 2')
xlabel('Avg Stride Dist')
ylabel('Avg Peak Mag')
hold off


end