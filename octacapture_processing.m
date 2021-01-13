%% processing two taps 12/2/19 Octa Capture - audacity
close all
filename = 'taps_120219.wav';
[y, Fs] = audioread(filename);
sig1 = y(:,1);
sig2 = y(:,2);
figure(1)
plot(sig1)
title('Left Sensor')
% time = linspace(0, length(sig1), length(sig1));
% figure(1)
% sig1 = sgolayfilt(sig1,1,71);% smoothdata(sig1, 'gaussian');
% % [pks, locs] = findpeaks(sig1);
% TF = islocalmax(sig1);
% % plot(time, sig1, time(locs), pks, 'or')
% plot(time,sig1,time(TF),sig1(TF),'r*')
figure(2)

plot(y(:,2))
title('Right Sensor')
figure(3)
% sig2 = sgolayfilt(sig2,1,71);
% plot(sig2)
plot(sig1-sig2)
figure(4)
plot(sig1+sig2)

%% 12/09/19 reading thru matlab 4 channels

close all

s = daq.createSession('directsound');

s.DurationInSeconds = 10.0;
% addAudioInputChannel(s,'Audio1','1');
% addAudioInputChannel(s,'Audio1','2');
% addAudioInputChannel(s,'Audio2','1');
% addAudioInputChannel(s,'Audio2','2');

addAudioInputChannel(s,'Audio5','1');
addAudioInputChannel(s,'Audio5','2');
addAudioInputChannel(s,'Audio6','1');
addAudioInputChannel(s,'Audio6','2');
data5 = startForeground(s);
%%
close all
figure(2)
plot(data5(:,1))
figure(3)
plot(data5(:,2))
figure(4)
plot(data5(:,3))
figure(5)
plot(data5(:,4))

%% 010620 identifying between 4 quadrants
close all
LT = data(:,3);
RT = data(:,1);
LB = data(:,4);
RB = data(:,2);

ntaps = 8;

LT_peaks = zeros(1,ntaps);
RT_peaks = zeros(1,ntaps);
LB_peaks = zeros(1,ntaps);
RB_peaks = zeros(1,ntaps);

for i = 1:ntaps
    % find time of index
    [~,idx] = max(abs(LT));
    LT_peaks(i) = idx;
    minrange = idx - 500;
    maxrange = idx + 12000;
    [~,idx] = max(abs(RT(minrange:maxrange)));
    RT_peaks(i) = idx+minrange;
    [~,idx] = max(abs(LB(minrange:maxrange)));
    LB_peaks(i) = idx+minrange;
    [~,idx] = max(abs(RB(minrange:maxrange)));
    RB_peaks(i) = idx+minrange;
    
    % eliminate this hit from array
    idx = idx + minrange;
    idx = idx - 500;
    idx = max(idx, 0);
    for j = 1:12000
        LT(idx+j) = 0;
        RT(idx+j) = 0;
        LB(idx+j) = 0;
        RB(idx+j) = 0;
    end
end

quadrants = zeros(1,ntaps);
for k = 1:ntaps
    [~, i] = max(LT_peaks); % max, so last hit to first
    right = RT_peaks(i) <= RB_peaks(i);
    top = RT_peaks(i) <= LT_peaks(i);
    left = LT_peaks(i) <= LB_peaks(i);
    bottom = RB_peaks(i) <= LB_peaks(i);
    if right && top && left && bottom
        hit = 1; % quadrant 1, counterclockwise from RT
    elseif ~right && top && ~left && bottom
        hit = 4;
    elseif right && ~top && left && ~bottom
        hit = 2;
    elseif ~right && ~top && ~left && ~bottom
        hit = 3;
    else
        hit = 0;
    end
    LT_peaks(i) = 0;
    quadrants(ntaps+1-k) = hit;
end

vidObj = VideoWriter('QuadrantsFind2_010620.mov');
open(vidObj);
figure(1)
xlim([-1 1])
ylim([-1 1])
hold on
for n = 1:ntaps
    xlim([-1 1])
    ylim([-1 1])
    quadrants(n)
    hold on
    x = linspace(-1,1,100);
    plot(0,x,'.-')
    plot(x,0,'.-')
    if quadrants(n) == 1
        plot(0.5,0.5,'*','MarkerSize', 25)
%         for x = linspace(0,1,30)
%             for y = linspace(0,1,30)
%                 plot(x,y,'*')
%             end
%         end
    elseif quadrants(n) == 2
        plot(-0.5,0.5,'*','MarkerSize', 25)
%         x = linspace(-1,0,30);
%         y = linspace(0,1,30);
%         plot(x,0,'*')
%         plot(x,1,'*')
%         plot(y,0,'*')
%         plot(y,-1,'*')
%         for x = linspace(-1,0,30)
%             for y = linspace(0,1,30)
%                 plot(x,y,'*')
%             end
%         end
    elseif quadrants(n) == 3
        plot(-0.5,-0.5,'*','MarkerSize', 25)
%         x = linspace(0,1,30);
%         y = linspace(0,1,30);
%         plot(x,0,'*')
%         plot(x,1,'*')
%         plot(y,0,'*')
%         plot(y,1,'*')
%         for x = linspace(-1,0,30)
%             for y = linspace(-1,0,30)
%                 plot(x,y,'*')
%             end
%         end
    elseif quadrants(n) == 4
        plot(0.5,-0.5,'*','MarkerSize', 25)
%         for x = linspace(0,1,30)
%             for y = linspace(-1,0,30)
%                 plot(x,y,'*')
%             end
%         end
    end
    hold off
    currFrame = getframe(gcf);
    for g = 1:30
        writeVideo(vidObj,currFrame);
    end
    clf
end
close(vidObj);

%% localization 021420
data = tips;
y = [];
x = [];
dLB = hilbert(data(:,LB));
dLT = hilbert(data(:,LT));
dRB = hilbert(data(:,RB));
dRT = hilbert(data(:,RT));
seconds = round(length(dLB)/1e5);
pkLB = [];
locLB = [];
pkLT = [];
locLT = [];
pkRB = [];
locRB = [];
pkRT = [];
locRT = [];

% AIC on 1 sec windows
for s = 0:(seconds-1)
    i_min = s*1e5+1;
    i_max = (s+1)*1e5;
    if s == seconds - 1
        i_max = length(dLB);
    end
    [pk,loc] = findpeaks(-1.*AIC_calc(dLB(i_min:i_max)));
    pkLB = pkLB + pk;
    locLB = locLB + loc;
    
    [pk,loc] = findpeaks(-1.*AIC_calc(dLT(i_min:i_max)));
    pkLT = pkLT + pk;
    locLT = locLT + loc;
    
    [pk,loc] = findpeaks(-1.*AIC_calc(dRB(i_min:i_max)));
    pkRB = pkRB + pk;
    locRB = locRB + loc;
    
    [pk,loc] = findpeaks(-1.*AIC_calc(dRT(i_min:i_max)));
    pkRT = pkRT + pk;
    locRT = locRT + loc;
%     [pk,loc] = findpeaks(-1.*AIC_calc(dLB(i_min:i_max)));
%     pkLB = pkLB + pk;
%     locLB = locLB + loc;
%     
%     [pk,loc] = findpeaks(-1.*AIC_calc(dLT(i_min:i_max)));
%     pkLT = pkLT + pk;
%     locLT = locLT + loc;
%     
%     [pk,loc] = findpeaks(-1.*AIC_calc(dRB(i_min:i_max)));
%     pkRB = pkRB + pk;
%     locRB = locRB + loc;
%     
%     [pk,loc] = findpeaks(-1.*AIC_calc(dRT(i_min:i_max)));
%     pkRT = pkRT + pk;
%     locRT = locRT + loc;
end

plot(AIC_calc(abs(dLT(i_min:i_max))) - smoothdata(AIC_calc(abs(dLT(i_min:i_max)))))
