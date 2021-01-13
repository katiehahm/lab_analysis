%% 6/10/20
% openExample('daq/demo_compactdaq_intro_acquisition')
clear all
% close all
% devices = daq.getDevices;
s = daq.createSession('ni');
s.addAnalogInputChannel('cDAQ1Mod1', 0, 'accelerometer');
s.addAnalogInputChannel('cDAQ1Mod1', 1, 'accelerometer');
s.addAnalogInputChannel('cDAQ1Mod1', 2, 'accelerometer');
% s.addAnalogInputChannel('cDAQ1Mod2', 0, 'accelerometer');
% s.addAnalogInputChannel('cDAQ1Mod2', 1, 'accelerometer');
% s.addAnalogInputChannel('cDAQ1Mod2', 2, 'accelerometer');
Fs = 12800;
disp('here')
s.Rate = Fs; % scan rate - scans/second
disp('here2')

%% 6/10/20
% ch2 = s.Channels(1);
ch1 = s.Channels(1);
ch2 = s.Channels(2);
ch3 = s.Channels(3);
% ch4 = s.Channels(4);
% ch5 = s.Channels(5);
% ch6 = s.Channels(6);

ch1.Sensitivity = 1; % V/g
ch2.Sensitivity = 1; % V/g
ch3.Sensitivity = 1; % V/g
% ch4.Sensitivity = 1; % V/g
% ch5.Sensitivity = 1; % V/g
% ch6.Sensitivity = 1; % V/g

ch1.ExcitationCurrent = 0.004; % 4 mA
ch2.ExcitationCurrent = 0.004; % 4 mA
ch3.ExcitationCurrent = 0.004; % 4 mA
% ch4.ExcitationCurrent = 0.004; % 4 mA
% ch5.ExcitationCurrent = 0.004; % 4 mA
% ch6.ExcitationCurrent = 0.004; % 4 mA

%% data plot testing
disp('Starting data collection')
s.DurationInSeconds = 90;
[datas, times] = s.startForeground;
disp('Data collection terminated');
filename = sprintf('data_%s', datestr(now,'mm-dd-yyyy_HH-MM'));
filename = append("data\", filename);
save(filename, 'datas', 'times')
disp(append("Saved as ", filename))

%% 6/16/20 plot continuous data
disp('Starting data collection')
s.DurationInSeconds = 1;
datas = [];
times = [];
[d, t] = s.startForeground;
datas = [datas; d];
times = [times, t];

disp('Press button to terminate')
figure;
ButtonHandle = uicontrol('Style', 'PushButton', ...
                         'String', 'Stop loop', ...
                         'Position', [30 25 500 370], ...
                         'Callback', 'delete(gcbf)');
% figure;
% hold on
for k = 1:1e6
    if ~ishandle(ButtonHandle)
        disp('Data collection terminated');
        break;
    else
        [d, t] = s.startForeground;
        datas = [datas; d];
        times = [times; t+k];
        disp(k)
%         if k > 100
%             figure;
%             hold on
%         end
%         hold on
%         subplot(3,1,1)
%         hold on
%         plot(t+k, d(:,1))
%         subplot(3,1,2)
%         hold on
%         plot(t+k, d(:,2))
%         subplot(3,1,3)
%         hold on
%         plot(t+k, d(:,3))
    end
end
%%
filename = sprintf('data_%s', datestr(now,'mm-dd-yyyy_HH-MM'));
filename = append("data\", filename);
save(filename, 'datas', 'times')
disp(append("Saved as ", filename))

% %%  6/16/20 plot 3 data
% figure;
% subplot(3,1,1)
% plot(times, datas(:,1))
% subplot(3,1,2)
% plot(times, datas(:,2))
% subplot(3,1,3)
% plot(times, datas(:,3))
% 
% xlabel('Time (secs)');
% ylabel('Voltage')
% 
% %% 6/10/20 plot all 6 data
% figure;
% subplot(3,2,1)
% plot(time, data(:,1))
% subplot(3,2,2)
% plot(time, data(:,2))
% subplot(3,2,3)
% plot(time, data(:,3))
% subplot(3,2,4)
% plot(time, data(:,4))
% subplot(3,2,5)
% plot(time, data(:,5))
% subplot(3,2,6)
% plot(time, data(:,6))
% % plot(time,data);
% xlabel('Time (secs)');
% ylabel('Voltage')