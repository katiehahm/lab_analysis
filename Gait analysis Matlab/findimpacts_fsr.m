function [impacts, Rheel, Lheel, Rtoe, Ltoe] = findimpacts_fsr(fsrTime,fsrData,Mfsr)
% finds the heel impacts from fsr data
% by looking at right and left heel strikes, findpeaks, and filtering out
% wide peaks and peaks with no =0 between them
% impacts = [heel peak idx, heel peak time, heel peak mag, 1(R) or 0(L)]
% R/Lheel = [initial impact index, init impact time, peak idx, peak magnitude]
% R/Ltoe = [end impact index, end impact time, peak idx, peak magnitude]
% 6/8/21, edited 9/1/21

Fs = 519;

figure;
plot(fsrData(:,Mfsr('Rheel')))
findpeaks(fsrData(:,Mfsr('Rheel')),'MinPeakProminence',10,'Annotate','extents','MinPeakDistance',0.5)
[pksR,locsR,widthsR,~] = findpeaks(fsrData(:,Mfsr('Rheel')),'MinPeakProminence',10,'Annotate','extents','MinPeakDistance',0.5); % 4th element is prominence
title('Right heel')

width_thresh = 500;
for i = 1:length(widthsR) % filter out wide peaks
    if widthsR(i) > width_thresh
        pksR(i) = 0;
        locsR(i) = 0;
    end
end
pksR = nonzeros(pksR);
locsR = nonzeros(locsR);

figure;
plot(fsrData(:,Mfsr('Lheel')))
findpeaks(fsrData(:,Mfsr('Lheel')),'MinPeakProminence',10,'Annotate','extents','MinPeakDistance',0.5)
[pksL,locsL,widthsL,~] = findpeaks(fsrData(:,Mfsr('Lheel')),'MinPeakProminence',10,'Annotate','extents','MinPeakDistance',0.5);
title('Left heel')

for i = 1:length(widthsL) % filter out wide peaks
    if widthsL(i) > width_thresh
        pksL(i) = 0;
        locsL(i) = 0;
    end
end
pksL = nonzeros(pksL);
locsL = nonzeros(locsL);

% filter out peaks that don't have heel go closer to 0
omitR = [];
for i = 2:length(pksR)
    btw_pk = fsrData(locsR(i-1):locsR(i),Mfsr('Rheel'));
    if ~any(btw_pk < 1)
        omitR(end+1) = i;
    end
end
locsR(omitR) = [];
pksR(omitR) = [];

omitL = [];
for i = 2:length(pksL)
    btw_pk = fsrData(locsL(i-1):locsL(i),Mfsr('Lheel'));
    if ~any(btw_pk < 1)
        omitL(end+1) = i;
    end
end
locsL(omitL) = [];
pksL(omitL) = [];

% find actual start of heel strike instead of peak
heel_to_peak = 0.2; % max 0.2s from onset of heel to peak
NL = length(locsL);
NR = length(locsR);
Rheel = zeros(NR, 4);
Lheel = zeros(NL, 4);
for i = 1:length(pksR)
    curr_idx = locsR(i);
    start_idx = curr_idx - round(heel_to_peak*Fs);
    heel_window = fsrData(start_idx:curr_idx,Mfsr('Rheel'));
    ipt = findchangepts(heel_window);
    Rheel(i,1) = ipt + start_idx - 1;
    Rheel(i,2) = fsrTime(Rheel(i,1));
end

Rheel(:,3) = locsR;
Rheel(:,4) = pksR;

for i = 1:length(pksL)
    curr_idx = locsL(i);
    start_idx = curr_idx - round(heel_to_peak*Fs);
    heel_window = fsrData(start_idx:curr_idx,Mfsr('Lheel'));
    ipt = findchangepts(heel_window);
    Lheel(i,1) = ipt + start_idx - 1;
    Lheel(i,2) = fsrTime(Lheel(i,1));
end

Lheel(:,3) = locsL;
Lheel(:,4) = pksL;

% sort them in time order
impacts = [locsR,fsrTime(locsR)',pksR,ones(length(locsR),1); locsL, fsrTime(locsL)',pksL, zeros(length(locsL),1)];
impacts = sortrows(impacts);

% find toe times
heel_to_toe_window = 1; % 1 sec window between heel peak and toe off
heel_to_toe_idx = heel_to_toe_window*Fs;
Ltoe = zeros(NL,4);
Rtoe = zeros(NR,4);
for i = 1:NL
    heel_idx = locsL(i);
    toe_window = fsrData(heel_idx:heel_idx + heel_to_toe_idx,Mfsr('Ltoe'));
    % toe peak
    [pk_mag, pk_idx] = max(toe_window);
    Ltoe(i,4) = pk_mag;
    Ltoe(i,3) = pk_idx + heel_idx - 1;
    % toe off
    ipt = findchangepts(toe_window);
    Ltoe(i,1) = ipt + heel_idx - 1;
    Ltoe(i,2) = fsrTime(Ltoe(i,1));
end

for i = 1:NR
    heel_idx = locsR(i);
    toe_window = fsrData(heel_idx:heel_idx+heel_to_toe_idx,Mfsr('Rtoe'));
    % toe peak
    [pk_mag, pk_idx] = max(toe_window);
    Rtoe(i,4) = pk_mag;
    Rtoe(i,3) = pk_idx + heel_idx - 1;
    % toe off
    ipt = findchangepts(toe_window);
    Rtoe(i,1) = ipt + heel_idx - 1;
    Rtoe(i,2) = fsrTime(Rtoe(i,1));
end

% visually check
figure; hold on
plot(fsrData(:,Mfsr('Rheel')))
plot(Rheel(:,1),0,'rx','MarkerSize',8)
plot(Rheel(:,3),Rheel(:,4),'bx','MarkerSize',8)
title('Right heel')

figure; hold on
plot(fsrData(:,Mfsr('Lheel')))
plot(Lheel(:,1),0,'rx','MarkerSize',8)
plot(Lheel(:,3),Lheel(:,4),'bx','MarkerSize',8)
title('Left heel')

figure; hold on
plot(fsrData(:,Mfsr('Ltoe')))
plot(Ltoe(:,1),0,'rx','MarkerSize',8)
plot(Ltoe(:,3),Ltoe(:,4),'bx','MarkerSize',8)
title('Left toe')

figure; hold on
plot(fsrData(:,Mfsr('Rtoe')))
plot(Rtoe(:,1),0,'rx','MarkerSize',8)
plot(Rtoe(:,3),Rtoe(:,4),'bx','MarkerSize',8)
title('Right toe')
end