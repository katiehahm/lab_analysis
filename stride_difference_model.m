load('data/data_09-25-2020_16-01.mat') %og model using 9-21 15-28
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

impactN = 13; % variable input value
[distanceval_table]=stridetimemodel(clean_data,impactN,Fs,loc_names)
%% 
function[distanceval_table]=stridetimemodel(clean_data,impactN,Fs,loc_names)
[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
sortedpeakidx=[sort(peak_idx(1,:));sort(peak_idx(2,:)); sort(peak_idx(3,:))]
peakdistdiff = [];
diff1=[]
diff2=[]
diff3=[]
for i = 2:impactN
    diff1 = [diff1 sortedpeakidx(1,i)-sortedpeakidx(1,i-1)];
    diff2 = [diff2 sortedpeakidx(2,i)-sortedpeakidx(2,i-1)];
    diff3 = [diff3 sortedpeakidx(3,i)-sortedpeakidx(3,i-1)];  
end
diff3
rowsToDelete = diff1 > 30000
diff1(rowsToDelete) = []
rowsToDelete = diff2 > 30000
diff2(rowsToDelete) = []
rowsToDelete = diff3 > 30000
diff3(rowsToDelete) = []
peakdistdiff = [diff1' diff2' diff3'];
 distance_foot= [peakdistdiff [1 2 1 2 1 2 1 2 1 2 1]'];
 distanceval_table = array2table(distance_foot,'VariableNames',{'Sensor1','Sensor2','Sensor3','foot'});


foot1 = [];
foot2=[];
for i = 1:11
    if distance_foot(i,4)== 1 & distance_foot(i,1)<30000
        foot1 = [foot1;distance_foot(i,1:3)];
    end
    if distance_foot(i,4)== 2 & distance_foot(i,1)<30000
        foot2 = [foot2;distance_foot(i,1:3)];
    end
end


x = categorical({'f1s1','f1s2','f1s3','f2s1', 'f2s2', 'f2s3'});
x = reordercats(x,{'f1s1','f1s2','f1s3','f2s1', 'f2s2', 'f2s3'});
y=[mean(foot1(:,1)) mean(foot1(:,2)) mean(foot1(:,3)) mean(foot2(:,1)) mean(foot2(:,2)) mean(foot2(:,3))];
figure
bar(x,y)
title('Average stride time difference from each sensor per foot')
ylabel('Stride Time (idx)')

diffx = categorical({'s1', 's2', 's3'});
diffx = reordercats(diffx,{'s1', 's2', 's3'});
diffy = [y(4)-y(1) y(5)-y(2) y(6)-y(3)];
figure
bar(diffx,diffy)
title('Foot 2 - Foot1')
ylabel('Stride Time (idx)')
end
