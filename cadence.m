%cluster and cadence per cluster
function [diffsinlength1,diffsinlength2,diffsinlength3] = cadence(peak_idx,impactN,Fs) 
sortedpeakidx= [sort(peak_idx(1,:));sort(peak_idx(2,:));sort(peak_idx(3,:))];
 diff1=[];
 diff2=[];
 diff3=[];
 for i = 2:impactN
    diff1 = [diff1 sortedpeakidx(1,i)-sortedpeakidx(1,i-1)];
    diff2 = [diff2 sortedpeakidx(2,i)-sortedpeakidx(2,i-1)];
    diff3 = [diff3 sortedpeakidx(3,i)-sortedpeakidx(3,i-1)];
 end
  diff1=diff1';
  diff2=diff2';
  diff3=diff3';
 diff1(diff1(:,1)>30000, :)= [];
 diff2(diff2(:,1)>30000, :)= [];
 diff3(diff3(:,1)>30000, :)= [];
cadence1 = (nnz(diff1)+1)/(sum(diff1(:,1))/Fs); %strikes per second
cadence2 = (nnz(diff1)+1)/(sum(diff2(:,1))/Fs);
cadence3 = (nnz(diff1)+1)/(sum(diff3(:,1))/Fs);
cadencetotal = [cadence1 cadence2 cadence3];
diffsinlength1=diff1/Fs;
diffsinlength2=diff2/Fs;
diffsinlength3=diff3/Fs;
% diff = [diff1 diff2 diff3];
% avgstridelengths = mean(diff)/Fs;
% avgstridelength = mean(avgstridelengths);
% stddevstridelength1= std(diff1)/Fs;
% stddevstridelength2= std(diff2)/Fs;
% stddevstridelength3= std(diff3)/Fs;
% stddevstridelength=[stddevstridelength1 stddevstridelength2 stddevstridelength3];
% figure
% bar(diff1/Fs)
% hold on
% plot(xlim,[avgstridelengths(1) avgstridelengths(1)], 'r')
% plot(xlim,[avgstridelengths(1)+stddevstridelength1 avgstridelengths(1)+stddevstridelength1], 'k')
% plot(xlim,[avgstridelengths(1)-stddevstridelength1 avgstridelengths(1)-stddevstridelength1], 'k')
% title('Sensor 1')
% ylim([0 1.8])
% hold off
% figure
% bar(diff2/Fs)
% hold on
% plot(xlim,[avgstridelengths(1,2) avgstridelengths(1,2)], 'r')
% plot(xlim,[avgstridelengths(1,2)+stddevstridelength2 avgstridelengths(2)+stddevstridelength2], 'k')
% plot(xlim,[avgstridelengths(1,2)-stddevstridelength2 avgstridelengths(2)-stddevstridelength2], 'k')
% title('Sensor 2')
% ylim([0 1.8])
% hold off
% figure
% bar(diff3/Fs)
% hold on
% plot(xlim,[avgstridelengths(1,3) avgstridelengths(1,3)], 'r')
% plot(xlim,[avgstridelengths(1,3)+stddevstridelength3 avgstridelengths(3)+stddevstridelength3], 'k')
% plot(xlim,[avgstridelengths(1,3)-stddevstridelength3 avgstridelengths(3)-stddevstridelength3], 'k')
% title('Sensor 3')
% ylim([0 1.8])
% hold off
end