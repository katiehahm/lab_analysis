function [impacts] = manual_fix_fsr(impacts,fsrData,Mfsr,heel_start_wrong,heel_start_right,heel_pk_wrong,heel_pk_right)
% takes in manual numbers to change for fsr heel strike and peak
% inputs are arrays
% 11/5/21

if ~isempty(heel_start_wrong)
    for i = 1:length(heel_start_wrong)
        indeces = abs(impacts(:,1) - heel_start_wrong(i));
        [~,change_idx] = min(indeces);
        impacts(change_idx,1) = heel_start_right(i);
    end
end

if ~isempty(heel_pk_wrong)
    for i = 1:length(heel_pk_wrong)
        indeces = abs(impacts(:,2) - heel_pk_wrong(i));
        [~,change_idx] = min(indeces);
        impacts(change_idx,2) = heel_pk_right(i);
    end
end

% visually check
RorL = length(impacts(1,:)); % always the last column
figure;
subplot(2,1,1)
hold on
plot(fsrData(:,Mfsr('Lheel')))
lefts = find(impacts(:,RorL) == 0);
plot(impacts(lefts,1),fsrData(impacts(lefts,1),Mfsr('Lheel')),'rx','MarkerSize',12)
plot(impacts(lefts,2),fsrData(impacts(lefts,2),Mfsr('Lheel')),'bx','MarkerSize',12)
title('Left heel')

subplot(2,1,2)
hold on
plot(fsrData(:,Mfsr('Rheel')))
rights = find(impacts(:,RorL) == 1);
plot(impacts(rights,1),fsrData(impacts(rights,1),Mfsr('Rheel')),'rx','MarkerSize',8)
plot(impacts(rights,2),fsrData(impacts(rights,2),Mfsr('Rheel')),'bx','MarkerSize',12)
title('Right heel')