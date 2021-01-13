function [peakdistdifftable] = stridetime_anaylsisfunction(peak_idx,impactN)

sortedpeakidx=[sort(peak_idx(1,:));sort(peak_idx(2,:)); sort(peak_idx(3,:))];
peakdistdiff = [];
diff1=[];
diff2=[];
diff3=[];
for i = 2:impactN
    diff1 = [diff1 sortedpeakidx(1,i)-sortedpeakidx(1,i-1)];
    diff2 = [diff2 sortedpeakidx(2,i)-sortedpeakidx(2,i-1)];
    diff3 = [diff3 sortedpeakidx(3,i)-sortedpeakidx(3,i-1)];  
end
peakdistdiff = [diff1' diff2' diff3' [1 2 1 2 1 2 1 2 1 2 1 2 1 2 1]'];
peakdistdifftable = array2table(peakdistdiff,'VariableNames',{'Sensor1','Sensor2','Sensor3','foot'});

% teststepsdata = DistancefootModel.predictFcn(peakdistdifftable);
% % teststepsdata=trainedClassifier.predictFcn(peakdistdifftable);
% peak_diffandstepinfo = [peakdistdiff teststepsdata];

% foot1 = [];
% foot2=[];
% for i = 1:impactN-1
%     if peak_diffandstepinfo(i,4)== 1 & peak_diffandstepinfo(i,1)<20000
%         foot1 = [foot1;peak_diffandstepinfo(i,1:3)];
%     end
%     if peak_diffandstepinfo(i,4)== 2 & peak_diffandstepinfo(i,1)<20000
%         foot2 = [foot2;peak_diffandstepinfo(i,1:3)];
%     end
end

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
%  end