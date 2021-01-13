load('data/data_09-21-2020_15-28.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

impactN = 24; % variable input value
clusterN = 6;


[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
% peakcomb = [peak_idx(:,1) peak_val];
% [~,idx] = sort(peakcomb(:,1)); % sort just the first column
% sortedmat = peakcomb(idx,:);  
% peak_val2 = sortedmat(:,2:4);
% foot_pattern = [peak_val [1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2]'];
% peakval_table = array2table(foot_pattern,'VariableNames',{'Sensor_1','Sensor_2','Sensor3','foot'});
% 
% foot1 = [];
% foot2=[];
% for i = 1:24
%     if foot_pattern(i,4)== 1
%         foot1 = [foot1;foot_pattern(i,1:3)];
%     else
%         foot2 = [foot2;foot_pattern(i,1:3)];
%     end
% end
% 
% 
% x = categorical({'f1s1','f1s2','f1s3','f2s1', 'f2s2', 'f2s3'});
% x = reordercats(x,{'f1s1','f1s2','f1s3','f2s1', 'f2s2', 'f2s3'});
% y=[mean(foot1(:,1)) mean(foot1(:,2)) mean(foot1(:,3)) mean(foot2(:,1)) mean(foot2(:,2)) mean(foot2(:,3))];
% figure
% bar(x,y)
% title('Average magnitude of peak from each sensor per foot')
% ylabel('Volts (V)')
% 
% diffx = categorical({'s1', 's2', 's3'});
% diffx = reordercats(diffx,{'s1', 's2', 's3'});
% diffy = [y(4)-y(1) y(5)-y(2) y(6)-y(3)];
% figure
% bar(diffx,diffy)
% title('Foot 2 - Foot1')
% ylabel('Volts (V)')

%  
% foot_clusterpattern = [peak_val2 [1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4]' [1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2]'];
% peakvalcluster_table = array2table(foot_clusterpattern,'VariableNames',{'Sensor_1','Sensor_2','Sensor3','cluster', 'foot'});

% avgpeak=[];
% for i = 1:24
%     avgpeak = [avgpeak mean(peak_val2(i,:))];
% end
% 
% 
% X=[avgpeak',avgpeak'];
% rng(1); % For reproducibility
% [idx,C]=kmeans(X,2);
% 
% figure;
% plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)
% hold on
% plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)
% plot(C(:,1),C(:,2),'kx',...
%      'MarkerSize',15,'LineWidth',3) 
% legend('Cluster 1','Cluster 2','Centroids',...
%        'Location','NW')
% title 'Cluster Assignments and Centroids'
% hold off