function [impacts] = manual_fix_fsr2pp(impacts,fsrData,heel_start_wrong,heel_start_right,heel_pk_wrong,heel_pk_right,delete_pks)
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

if ~isempty(delete_pks)
    for i = 1:length(delete_pks)
        indeces = abs(impacts(:,2) - delete_pks(i));
        [~,change_idx] = min(indeces);
        impacts(change_idx,:) = [];
    end
end

% visually check
for i = 1:4
    figure;
    idxs = find(impacts(:,4) == i);
    hold on
    plot(fsrData(:,i))
    plot(impacts(idxs,1),fsrData(impacts(idxs,1),i), 'rx', 'MarkerSize', 12)
    plot(impacts(idxs,2),impacts(idxs,3),'bx','MarkerSize',12)
    title(['Foot ', string(i)])
end