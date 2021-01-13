%% testing function
load('data/data_09-25-2020_16-01.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);
[testedmodel]=foottest(trainedClassifier,clean_data,impactN,Fs,loc_names)

%% training data (must be run first to generate classifier variable)

load('data/data_09-25-2020_16-01.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

impactN = 13; % variable input value
[peakval_table trainedClassifier]=footmode(clean_data,impactN,Fs,loc_names)

function [peakvalcluster_test trainedClassifier]=footmode(clean_data,impactN,Fs,loc_names)
    [onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
    peak1 = [peak_idx(1,:); peak_val(1,:)];
    peak2 = [peak_idx(2,:); peak_val(2,:)];
    peak3 = [peak_idx(3,:); peak_val(3,:)];
    peak1=peak1';
    peak2=peak2';
    peak3=peak3';
    [~,idx] = sort(peak1(:,1)); % sort just the first column
    sortedpeak1 = peak1(idx,:);
    [~,idx] = sort(peak2(:,1)); % sort just the first column
    sortedpeak2 = peak2(idx,:);
    [~,idx] = sort(peak3(:,1)); % sort just the first column
    sortedpeak3 = peak3(idx,:);
    % foot_pattern = [peak_val [1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2]'];
    % peakval_table = array2table(foot_pattern,'VariableNames',{'Sensor_1','Sensor_2','Sensor3','foot'});
    %
    
    foot_clusterpatterntest = [sortedpeak1(:,2) sortedpeak2(:,2) sortedpeak3(:,2) [1 2 1 2 1 2 1 2 1 2 1 2 1]'];
    peakvalcluster_test = array2table(foot_clusterpatterntest,'VariableNames',{'Sensor_1','Sensor_2','Sensor3','foot'});
    [trainedClassifier, validationAccuracy] = trainClassifierfoot(peakvalcluster_test)
%     foot1 = [];
%     foot2=[];
%     for i = 1:im
%         if peak_magandstepinfo(i,4)== 1
%             foot1 = [foot1;peak_magandstepinfo(i,1:3)];
%         else
%             foot2 = [foot2;peak_magandstepinfo(i,1:3)];
%         end
%     end
%     
%     
%     x = categorical({'f1s1','f1s2','f1s3','f2s1', 'f2s2', 'f2s3'});
%     x = reordercats(x,{'f1s1','f1s2','f1s3','f2s1', 'f2s2', 'f2s3'});
%     y=[mean(foot1(:,1)) mean(foot1(:,2)) mean(foot1(:,3)) mean(foot2(:,1)) mean(foot2(:,2)) mean(foot2(:,3))];
%     figure
%     bar(x,y)
%     title('Average magnitude of peak from each sensor per foot')
%     ylabel('Volts (V)')
%     
%     diffx = categorical({'s1', 's2', 's3'});
%     diffx = reordercats(diffx,{'s1', 's2', 's3'});
%     diffy = [y(4)-y(1) y(5)-y(2) y(6)-y(3)];
%     figure
%     bar(diffx,diffy)
%     title('Foot 2 - Foot1')
%     ylabel('Volts (V)')
%     
%     % testdata1601 = trainedModelcluster.predictFcn(peakvalcluster_test);
%     
%     % peak_magandstepinfo = [peak_val2 testdata1601];
%     
%     foot1 = [];
%     foot2=[];
%     for i = 1:13
%         if peak_magandstepinfo(i,4)== 1
%             foot1 = [foot1;peak_magandstepinfo(i,1:3)];
%         else
%             foot2 = [foot2;peak_magandstepinfo(i,1:3)];
%         end
%     end
%     
%     
%     x = categorical({'f1s1','f1s2','f1s3','f2s1', 'f2s2', 'f2s3'});
%     x = reordercats(x,{'f1s1','f1s2','f1s3','f2s1', 'f2s2', 'f2s3'});
%     y=[mean(foot1(:,1)) mean(foot1(:,2)) mean(foot1(:,3)) mean(foot2(:,1)) mean(foot2(:,2)) mean(foot2(:,3))];
%     figure
%     bar(x,y)
%     title('Average magnitude of peak from each sensor per foot')
%     ylabel('Volts (V)')
%     
%     diffx = categorical({'s1', 's2', 's3'});
%     diffx = reordercats(diffx,{'s1', 's2', 's3'});
%     diffy = [y(4)-y(1) y(5)-y(2) y(6)-y(3)];
%     figure
%     bar(diffx,diffy)
%     title('Foot 2 - Foot1')
%     ylabel('Volts (V)')
end
        