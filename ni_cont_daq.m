%% 1/8/2021
clear all
close all
global datas
global times
Fs = 12800;

s = daq.createSession('ni');

% ch(1) = addAnalogInputChannel(s,'cDAQ1Mod1', 0, 'accelerometer');
% ch(2) = addAnalogInputChannel(s,'cDAQ1Mod1', 1, 'accelerometer');
% ch(3) = addAnalogInputChannel(s,'cDAQ1Mod1', 2, 'accelerometer');
% ch(4) = addAnalogInputChannel(s,'cDAQ1Mod2', 0, 'accelerometer');

ch1 = addAnalogInputChannel(s,'cDAQ1Mod1', 0, 'accelerometer');
ch2 = addAnalogInputChannel(s,'cDAQ1Mod1', 1, 'accelerometer');
ch3 = addAnalogInputChannel(s,'cDAQ1Mod1', 2, 'accelerometer');
ch4 = addAnalogInputChannel(s,'cDAQ1Mod2', 0, 'accelerometer');

ch1.Sensitivity = 1; % V/g 
ch2.Sensitivity = 1; % V/g
ch3.Sensitivity = 1; % V/g
ch4.Sensitivity = 1; % V/g

% ch.Sensitivity = 1; % V/g
lh = addlistener(s,'DataAvailable', @continuousAcq);

s.Rate = Fs;
s.NotifyWhenDataAvailableExceeds = Fs/4; %call ScansAvailableFcn twice/sec
s.IsContinuous = true;

% s.ScansAvailableFcn = @(src,event) closePlot(src,event);
% s.ScansAvailableFcnCount = Fs*10; % 10 secs

s.startBackground()
disp('Type stop(s) to terminate')
figure('Renderer', 'painters', 'Position', [900 10 1000 800])
hold on

% while s.IsRunning
%     pause(10)
%     close all
% end

%%

filename = sprintf('data_%s', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
filename = append("data\", filename);
save(filename, 'datas', 'times')
disp(append("Saved as ", filename))

%% using DataAcquisition; doesn't work yet 1/10/21
% https://www.mathworks.com/help/daq/daq.interfaces.dataacquisition.html
% https://www.mathworks.com/help/daq/transition-your-code-from-session-to-dataacquisition-interface.html

clear all
close all
global datas
global times
Fs = 12800;

dq = daq('ni');
addinput(dq, 'cDAQ1Mod1', 0, 'accelerometer');
addinput(dq, 'cDAQ1Mod1', 1, 'accelerometer');
addinput(dq, 'cDAQ1Mod1', 2, 'accelerometer');
addinput(dq, 'cDAQ1Mod2', 0, 'accelerometer');

dq.Rate = Fs;

s.ScansAvailableFcn = @(src,event) continuousAcq(src,event);
s.ScansAvailableFcnCount = Fs/2; % 10 secs

while dq.Running
    pause(10)
    close all
end

start(dq, "continuous");


% % ch(1) = addAnalogInputChannel(s,'cDAQ1Mod1', 0, 'accelerometer');
% % ch(2) = addAnalogInputChannel(s,'cDAQ1Mod1', 1, 'accelerometer');
% % ch(3) = addAnalogInputChannel(s,'cDAQ1Mod1', 2, 'accelerometer');
% % ch(4) = addAnalogInputChannel(s,'cDAQ1Mod2', 0, 'accelerometer');
% 
% ch1 = addAnalogInputChannel(s,'cDAQ1Mod1', 0, 'accelerometer');
% ch2 = addAnalogInputChannel(s,'cDAQ1Mod1', 1, 'accelerometer');
% ch3 = addAnalogInputChannel(s,'cDAQ1Mod1', 2, 'accelerometer');
% ch4 = addAnalogInputChannel(s,'cDAQ1Mod2', 0, 'accelerometer');
% 
% ch1.Sensitivity = 1; % V/g 
% ch2.Sensitivity = 1; % V/g
% ch3.Sensitivity = 1; % V/g
% ch4.Sensitivity = 1; % V/g
% 
% % ch.Sensitivity = 1; % V/g
% lh = addlistener(s,'DataAvailable', @continuousAcq);
% 
% s.Rate = Fs;
% s.NotifyWhenDataAvailableExceeds = Fs/2; %call ScansAvailableFcn twice/sec
% s.IsContinuous = true;
% 
% s.ScansAvailableFcn = @(src,event) closePlot(src,event);
% s.ScansAvailableFcnCount = Fs*10; % 10 secs
% 
% s.startBackground()
disp('Type stop(s) to terminate')
figure('Renderer', 'painters', 'Position', [900 10 1000 800])
hold on

% while s3.IsRunning
%     if 
% end
