function [impactT, impactsM, RheelT, RtoeT, LheelT, LtoeT] = findimpacts_fsr(fsrT,fsrD,Mfsr)
% finds the heel impacts from fsr data
% by looking at right and left heel strikes, findpeaks, and filtering out
% wide peaks and peaks with no =0 between them
% 6/8/21

figure;
plot(fsrD(:,Mfsr('Rheel')))
title('Right heel')
findpeaks(fsrD(:,Mfsr('Rheel')),'MinPeakProminence',10,'Annotate','extents','MinPeakDistance',0.5)
[pksR,locsR,widthsR,~] = findpeaks(fsrD(:,Mfsr('Rheel')),'MinPeakProminence',10,'Annotate','extents','MinPeakDistance',0.5); % 4th element is prominence

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
plot(fsrD(:,Mfsr('Lheel')))
title('Left heel')
findpeaks(fsrD(:,Mfsr('Lheel')),'MinPeakProminence',10,'Annotate','extents','MinPeakDistance',0.5)
[pksL,locsL,widthsL,~] = findpeaks(fsrD(:,Mfsr('Lheel')),'MinPeakProminence',10,'Annotate','extents','MinPeakDistance',0.5);

for i = 1:length(widthsL) % filter out wide peaks
    if widthsL(i) > width_thresh
        pksL(i) = 0;
        locsL(i) = 0;
    end
end
pksL = nonzeros(pksL);
locsL = nonzeros(locsL);
RheelT = [];
RtoeT = [];

impacts = [locsR(1),pksR(1);
        locsL(1),pksL(1)];
for i = 2:length(pksR)
    btw_pk = fsrD(locsR(i-1):locsR(i),Mfsr('Rheel'));
    if any(btw_pk < 1)
        impacts(end+1,:) = [locsR(i),pksR(i)];
    end
end
for i = 2:length(pksL)
    btw_pk = fsrD(locsL(i-1):locsL(i),Mfsr('Lheel'));
    if any(btw_pk < 1)
        impacts(end+1,:) = [locsL(i),pksL(i)];
    end
end
impacts = sortrows(impacts);

% visually check
figure; hold on
plot(fsrD(:,Mfsr('Lheel')))
plot(fsrD(:,Mfsr('Rheel')))
plot(impacts(:,1),impacts(:,2), 'ko','MarkerSize',10)

impactT = fsrT(impacts(:,1));

end

