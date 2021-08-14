clear all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\';

% ################# change ########################
% datestr = '05-09-2021_15-08-29';
% datestr = '06-04-2021_13-47-54';
% datestr = '06-10-2021_13-23-05';
% datestr = '06-10-2021_13-38-16';
% datestr = '06-10-2021_14-00-15';
datestr = '06-10-2021_14-17-29';
load([data_root_katie, 'ProcessedData\', datestr]) % ### change
% keySet = {'Rballi','Rtoe','Rheel','Rballo','Ltoe','Lheel','Lballi','Lballo'};
figure;
plot(fsrT,fsrD(:,3))
figure;
plot(fsrT,fsrD(:,6))