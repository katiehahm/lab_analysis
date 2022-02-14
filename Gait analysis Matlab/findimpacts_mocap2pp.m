function [coordinates] = findimpacts_mocap2pp(impacts,fsrTime,mocapT,allmocap,visualize)
% finds the locations from mocap 
% finds where the y-loc dips
% then confirms based on impacts found through fsr
% saves locations to coordinates matrix
% whichfoot: right is 1, left is 0
% 9/9/21, modified from before when it just used fsr to find the next coord
% now it inverts the y-sig, find peak, plus buffer, then store coord if ok
% with fsr
% extracted_idx_R/L stores the index of the extracted coordinate

% impactT = fsrTime(impacts(:,1));
coordinates = zeros(length(impacts),5); % x, z, person id
window = 80; % 80 samples of mocap to look for peak
leg_still = 20; % leg stands still for at least 20 pts

for i = 1:length(impacts)
    personID = impacts(i,4);
    mocap_avg = movmean(-1.*allmocap(:,2,personID), 10);
    heelTime = fsrTime(impacts(i,1));
    mocap_i = findTindex(heelTime, mocapT);
    clip = mocap_avg(mocap_i - window*0.5:mocap_i + window*0.5);
    [~,maxidx] = max(clip);
    coord_i = maxidx + mocap_i - window*0.5 + leg_still;
    coordinates(i,:) = [allmocap(coord_i,1,personID), allmocap(coord_i,2,personID), allmocap(coord_i,3,personID), coord_i, personID];
end

if visualize
    xyz = ["x", "y", "z"];
    for n = 1:4
        figure;
        person_idx = find(coordinates(:,5) == n);
        for k = 1:3
            subplot(3,1,k)
            hold on
            plot(allmocap(:,k,n))
            hold on
            plot(coordinates(person_idx,4), coordinates(person_idx,k), 'r.', 'MarkerSize', 10)
            title(['Foot ', string(n), ' coordinate ', xyz(k)])
        end
    end
end


end

