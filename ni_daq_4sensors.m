%% 12/10/20 acquire background data
% **** Run Section **** to save data to .mat
clear all
close all
global datas
global times
Fs = 12800;
s = daq.createSession('ni');
addAnalogInputChannel(s,'cDAQ1Mod1', 0, 'accelerometer');
addAnalogInputChannel(s,'cDAQ1Mod1', 1, 'accelerometer');
addAnalogInputChannel(s,'cDAQ1Mod1', 2, 'accelerometer');
addAnalogInputChannel(s,'cDAQ1Mod2', 0, 'accelerometer');

ch1 = s.Channels(1);
ch2 = s.Channels(2);
ch3 = s.Channels(3);
ch4 = s.Channels(4);
ch1.Sensitivity = 1; % V/g 
ch2.Sensitivity = 1; % V/g
ch3.Sensitivity = 1; % V/g
ch4.Sensitivity = 1; % V/g

s.Rate = Fs;
lh = addlistener(s,'DataAvailable', @continuousAcq);

s.NotifyWhenDataAvailableExceeds = Fs/2; %call ScansAvailableFcn twice/sec
s.IsContinuous = true;
s.startBackground()
disp('Type stop(s) to terminate')
figure('Renderer', 'painters', 'Position', [900 10 1000 800])
hold on
%% 7/8/20 for saving data
filename = sprintf('data_%s', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
filename = append("data\", filename);
save(filename, 'datas', 'times')
disp(append("Saved as ", filename))

