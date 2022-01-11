function [impacts] = findimpacts_fsr_accel(fsrTime,fsrData,Mfsr,L_dist,R_dist,min_threshL,min_threshR)
% finds the heel impacts from fsr data that includes accel
% by looking at right and left heel strikes, findpeaks, and filtering out
% wide peaks and peaks with no =0 between them
% impacts = [heel start idx, heel pk idx, heel pk mag, 1(R) or 0(L)]
% 9/15/21

[pksR,locsR,~,~] = findpeaks(fsrData(:,Mfsr('Rheel')),'MinPeakProminence',8,'Annotate','extents','MinPeakDistance',R_dist); % 4th element is prominence
[pksL,locsL,~,~] = findpeaks(fsrData(:,Mfsr('Lheel')),'MinPeakProminence',9,'Annotate','extents','MinPeakDistance',L_dist);

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

% always throw out first peak because often not clean
Rheel(1,:) = [];
Lheel(1,:) = [];
% sort them in time order
impacts = [Rheel(:,1),locsR(2:end),pksR(2:end),ones(length(locsR)-1,1);...
    Lheel(:,1), locsL(2:end), pksL(2:end), zeros(length(locsL)-1,1)];
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
end