%% Trying to find peaks without # impacts
data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_09-21-2020_15-27.mat').datas;  % This is one set of step data
data_1 = data(:,1);                                                                         % Only data from one sensor
[peaks, locs, w, p] = findpeaks(data(:,1), 'MinPeakHeight', 0.011, 'MinPeakDistance', 2000);

%% Testing process on more steps 
data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_09-21-2020_15-28.mat').datas;  % This is one set of step data
data_1 = data(:,1);                                                                         % Only data from one sensor
[peaks, locs, w, p] = findpeaks(data(:,1), 'MinPeakHeight', 0.011, 'MinPeakDistance', 2000);

%% Last test on step finder (using stomps this time)
data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_09-25-2020_15-57.mat').datas;  % This is one set of step data
data_1 = data(:,1);                                                                         % Only data from one sensor
[peaks, locs, w, p] = findpeaks(data(:,1), 'MinPeakHeight', 0.011, 'MinPeakDistance', 2000);

%% Testing TDOA3 function
data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_09-25-2020_15-57.mat').datas;
sensor1 = data(:,1);
[peaks, locs, w, p] = TDOA3(sensor1);