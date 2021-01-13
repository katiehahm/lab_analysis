%% Identifying peaks for 4 sensors (1-3)
data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_12-10-2020_12-16-21.mat').datas; % black ball 31.5
sens_1 = data(:, 1);
sens_2 = data(:, 2);
sens_3 = data(:, 3);
sens_4 = data(:, 4);
[peaks_1, idx_1, width_1, prominence_1] = TDOA3(sens_1);
[peaks_2, idx_2, width_2, prominence_2] = TDOA3(sens_2);
[peaks_3, idx_3, width_3, prominence_3] = TDOA3(sens_3);
[peaks_4, idx_4, width_4, prominence_4] = TDOA3(sens_4);
% Pos 1
avg_pos_1_sens_1 = mean(peaks_1(1:5,:),1);
avg_pos_1_sens_2 = mean(peaks_2(1:5,:),1);
avg_pos_1_sens_3 = mean(peaks_3(1:5,:),1);
avg_pos_1_sens_4 = mean(peaks_4(1:5,:),1);
% Pos 2
avg_pos_2_sens_1 = mean(peaks_1(6:10,:),1);
avg_pos_2_sens_2 = mean(peaks_2(6:10,:),1);
avg_pos_2_sens_3 = mean(peaks_3(6:10,:),1);
avg_pos_2_sens_4 = mean(peaks_4(6:10,:),1);
% Pos 3
avg_pos_3_sens_1 = mean(peaks_1(11:15,:),1);
avg_pos_3_sens_2 = mean(peaks_2(11:15,:),1);
avg_pos_3_sens_3 = mean(peaks_3(11:15,:),1);
avg_pos_3_sens_4 = mean(peaks_4(11:15,:),1);

%% Identifying peaks for 4 sensors (4-6)
data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_12-10-2020_12-18-11.mat').datas; % black ball 31.5
sens_1 = data(:, 1);
sens_2 = data(:, 2);
sens_3 = data(:, 3);
sens_4 = data(:, 4);
[peaks_1, idx_1, width_1, prominence_1] = TDOA3(sens_1);
[peaks_2, idx_2, width_2, prominence_2] = TDOA3(sens_2);
[peaks_3, idx_3, width_3, prominence_3] = TDOA3(sens_3);
[peaks_4, idx_4, width_4, prominence_4] = TDOA3(sens_4);
% Pos 4
avg_pos_4_sens_1 = mean(peaks_1(1:5,:),1);
avg_pos_4_sens_2 = mean(peaks_2(1:5,:),1);
avg_pos_4_sens_3 = mean(peaks_3(1:5,:),1);
avg_pos_4_sens_4 = mean(peaks_4(1:5,:),1);
% Pos 5
avg_pos_5_sens_1 = mean(peaks_1(6:10,:),1);
avg_pos_5_sens_2 = mean(peaks_2(6:10,:),1);
avg_pos_5_sens_3 = mean(peaks_3(6:10,:),1);
avg_pos_5_sens_4 = mean(peaks_4(6:10,:),1);
% Pos 6
avg_pos_6_sens_1 = mean(peaks_1(11:15,:),1);
avg_pos_6_sens_2 = mean(peaks_2(11:15,:),1);
avg_pos_6_sens_3 = mean(peaks_3(11:15,:),1);
avg_pos_6_sens_4 = mean(peaks_4(11:15,:),1);



%% Distances to sensors
% Pos 1
%   125" from sensor 1
%   65" from sensor 2
%   77" from sensor 3
%   155" from sensor 4

% Pos 2
%   115" from sensor 1
%   77" from sensor 2
%   87" from sensor 3
%   122" from sensor 4

% Pos 3
%   98" from sensor 1
%   91" from sensor 2
%   100" from sensor 3
%   106" from sensor 4

% Pos 4
%   83" from sensor 1
%   107" from sensor 2
%   115" from sensor 3
%   92" from sensor 4

% Pos 5
%   70" from sensor 1
%   124" from sensor 2
%   131" from sensor 3
%   81" from sensor 4

% Pos 6
%   61" from sensor 1
%   142" from sensor 2
%   148" from sensor 3
%   73" from sensor 4

%% Results from position 1


