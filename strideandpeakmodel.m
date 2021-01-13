load('data/data_10-02-2020_19-06.mat') %og model using 9-21 15-28
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

impactN = 20; % variable input value


[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
peaktotal=[peak_idx' peak_val'];
[~,idx] = sort(peaktotal(:,1)); % sort just the first column
sortedpeaktotal = peaktotal(idx,:);   % sort the whole matrix using the sort indices
%sortedpeakidx=[sort(peak_idx(1,:));sort(peak_idx(2,:)); sort(peak_idx(3,:))];

 peakdistdiff = [];
 diff1=[];
 diff2=[];
 diff3=[];
 for i = 2:impactN
    diff1 = [diff1 sortedpeaktotal(i,1)-sortedpeaktotal(i-1,1)];
    diff2 = [diff2 sortedpeaktotal(i,2)-sortedpeaktotal(i-1,2)];
    diff3 = [diff3 sortedpeaktotal(i,3)-sortedpeaktotal(i-1,3)];  
 end
 diffpeak1=[diff1' sortedpeaktotal(2:impactN,4)];
 diffpeak2=[diff2' sortedpeaktotal(2:impactN,5)];
 diffpeak3=[diff1' sortedpeaktotal(2:impactN,6)];
 j=impactN-1;
diffpeak1(diffpeak1(:,1)>20000, :)= [];
diffpeak2(diffpeak2(:,1)>20000, :)= [];
diffpeak3(diffpeak3(:,1)>20000, :)= [];

peakdistdiff = [diffpeak1 diffpeak2 diffpeak3];
distance_foot= [peakdistdiff [1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2]'];
distanceval_table = array2table(distance_foot,'VariableNames',{'Sensor1idx','Sensor1peak','Sensor2idx','Sensor2peak','Sensor3idx','Sensor3peak','foot'});
%  
% 
% foot1 = [];
% foot2=[];
% for i = 1:12
%     if distance_foot(i,4)== 1 & distance_foot(i,1)<20000
%         foot1 = [foot1;distance_foot(i,1:3)];
%     end
%     if distance_foot(i,4)== 2 & distance_foot(i,1)<20000
%         foot2 = [foot2;distance_foot(i,1:3)];
%     end
% end
% 
% 
% x = categorical({'f1s1','f1s2','f1s3','f2s1', 'f2s2', 'f2s3'});
% x = reordercats(x,{'f1s1','f1s2','f1s3','f2s1', 'f2s2', 'f2s3'});
% y=[mean(foot1(:,1)) mean(foot1(:,2)) mean(foot1(:,3)) mean(foot2(:,1)) mean(foot2(:,2)) mean(foot2(:,3))];
% figure
% bar(x,y)
% title('Average stride time difference from each sensor per foot')
% ylabel('Stride Time (idx)')
% 
% diffx = categorical({'s1', 's2', 's3'});
% diffx = reordercats(diffx,{'s1', 's2', 's3'});
% diffy = [y(4)-y(1) y(5)-y(2) y(6)-y(3)];
% figure
% bar(diffx,diffy)
% title('Foot 2 - Foot1')
% ylabel('Stride Time (idx)')
