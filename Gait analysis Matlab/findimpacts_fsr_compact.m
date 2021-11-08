function [impacts] = findimpacts_fsr_compact(fsrTime,fsrData,Mfsr,L_dist,R_dist,min_threshL,min_threshR,toe_threshL,toe_threshR)
% finds the heel & toe impacts from fsr data
% by looking at right and left heel strikes, findpeaks, and filtering out
% wide peaks and peaks with no =0 between them
% impacts = [heel start idx, heel pk idx, heel pk mag, toe end idx, toe pk idx, toe pk mag, 1(R) or 0(L)]
% 9/15/21

Fs = 518.5;

% figure;
% plot(fsrData(:,Mfsr('Rheel')))
[pksR,locsR,~,~] = findpeaks(fsrData(:,Mfsr('Rheel')),'MinPeakProminence',9,'Annotate','extents','MinPeakDistance',R_dist); % 4th element is prominence
% hold on
% plot(locsR,pksR,'bx')
% title('Right heel')

% figure;
% plot(fsrData(:,Mfsr('Lheel')))
[pksL,locsL,~,~] = findpeaks(fsrData(:,Mfsr('Lheel')),'MinPeakProminence',9,'Annotate','extents','MinPeakDistance',L_dist);
% hold on
% plot(locsL,pksL,'bx')
% title('Left heel')

% just for subj 4 count:
% pksL = pksL(4:end);
% locsL = locsL(4:end);
% widthsL = widthsL(4:end);

% filter out peaks that don't have heel go close to 0 in between
omitR = [];
for i = 2:length(pksR)
    btw_pk = fsrData(locsR(i-1):locsR(i),Mfsr('Rheel'));
    if ~any(btw_pk < min_threshR)
        omitR(end+1) = i;
    end
end
locsR(omitR) = [];
pksR(omitR) = [];

omitL = [];
for i = 2:length(pksL)
    btw_pk = fsrData(locsL(i-1):locsL(i),Mfsr('Lheel'));
    if ~any(btw_pk < min_threshL)
        omitL(end+1) = i;
    end
end
locsL(omitL) = [];
pksL(omitL) = [];

% find actual start of heel strike instead of peak
NL = length(locsL);
NR = length(locsR);
Rheel = zeros(NR, 4);
Lheel = zeros(NL, 4);
for i = 2:length(pksR)
    window = fsrData(locsR(i-1):locsR(i),Mfsr('Rheel'));
    [~,idx] = min(window);
    window = movmean(window,20);
%     dgrad = gradient(gradient(window));
    start_idx = aic_pick(window(idx:end),'to_peak') + idx + locsR(i-1);
    Rheel(i,1) = start_idx;
    Rheel(i,2) = fsrTime(start_idx);
end
Rheel(:,3) = locsR;
Rheel(:,4) = pksR;

for i = 2:length(pksL)
    window = fsrData(locsL(i-1):locsL(i),Mfsr('Lheel'));
    [~,idx] = min(window);
    window = movmean(window,20);
%     dgrad = gradient(gradient(window));
    start_idx = aic_pick(window(idx:end),'to_peak') + idx + locsL(i-1);
    Lheel(i,1) = start_idx;
    Lheel(i,2) = fsrTime(start_idx);
end
Lheel(:,3) = locsL;
Lheel(:,4) = pksL;

% find toe times
Ltoe = zeros(NL,4);
Rtoe = zeros(NR,4);
for i = 2:NL
    heel_idx = Lheel(i,1);
    if i == NL
        end_window = length(fsrData(:,1));
    else
        end_window = Lheel(i+1,1);
    end
    toe_window = fsrData(heel_idx:end_window,Mfsr('Ltoe'));
    % toe peak
    [pk_mag, pk_idx] = max(toe_window);
    Ltoe(i,4) = pk_mag;
    Ltoe(i,3) = pk_idx + heel_idx - 1;
    % toe off
    indeces = find(toe_window(pk_idx:end) < toe_threshL);
    if length(indeces) < 1 % couldn't find toe off, something wrong
        heel_idx
        figure;
        plot(toe_window(pk_idx:end))
    else
    end_idx = indeces(1) + pk_idx + heel_idx - 1;
    Ltoe(i,1) = end_idx;
    Ltoe(i,2) = fsrTime(end_idx);
    end
end

for i = 2:NR
    heel_idx = Rheel(i,1);
    if i == NR
        end_window = length(fsrData(:,1));
    else
        end_window = Rheel(i+1,1);
    end
    toe_window = fsrData(heel_idx:end_window,Mfsr('Rtoe'));
    % toe peak
    [pk_mag, pk_idx] = max(toe_window);
    Rtoe(i,4) = pk_mag;
    Rtoe(i,3) = pk_idx + heel_idx - 1;
    % toe off
    indeces = find(toe_window(pk_idx:end) < toe_threshR);
    if length(indeces) < 1 % couldn't find toe off, something wrong
        heel_idx
        figure;
        plot(toe_window(pk_idx:end))
    else
    end_idx = indeces(1) + pk_idx + heel_idx - 1;
    Rtoe(i,1) = end_idx;
    Rtoe(i,2) = fsrTime(end_idx);
    end
end

% always throw out first peak because often not clean
Rheel(1,:) = [];
Rtoe(1,:) = [];
Lheel(1,:) = [];
Ltoe(1,:) = [];
% sort them in time order
impacts = [Rheel(:,1),locsR(2:end),pksR(2:end),Rtoe(:,1),Rtoe(:,3),Rtoe(:,4),ones(length(locsR)-1,1);...
    Lheel(:,1), locsL(2:end), pksL(2:end), Ltoe(:,1),Ltoe(:,3),Ltoe(:,4),zeros(length(locsL)-1,1)];
impacts = sortrows(impacts,1);

% visually check
figure;
subplot(2,1,1)
hold on
plot(fsrData(:,Mfsr('Lheel')))
plot(Lheel(:,1),fsrData(Lheel(:,1),Mfsr('Lheel')),'rx','MarkerSize',12)
plot(Lheel(:,3),Lheel(:,4),'bx','MarkerSize',12)
title('Left heel')

subplot(2,1,2)
hold on
plot(fsrData(:,Mfsr('Rheel')))
plot(Rheel(:,1),fsrData(Rheel(:,1),Mfsr('Rheel')),'rx','MarkerSize',12)
plot(Rheel(:,3),Rheel(:,4),'bx','MarkerSize',12)
title('Right heel')

figure;
subplot(2,1,1)
hold on
plot(fsrData(:,Mfsr('Ltoe')))
plot(Ltoe(:,1),0,'rx','MarkerSize',12)
plot(Ltoe(:,3),Ltoe(:,4),'bx','MarkerSize',12)
title('Left toe')

subplot(2,1,2)
hold on
plot(fsrData(:,Mfsr('Rtoe')))
plot(Rtoe(:,1),0,'rx','MarkerSize',12)
plot(Rtoe(:,3),Rtoe(:,4),'bx','MarkerSize',12)
title('Right toe')
end