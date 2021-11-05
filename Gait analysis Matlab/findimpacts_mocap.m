function [extracted_idx_R, extracted_idx_L, coordinates, whichfoot] = findimpacts_mocap(impacts,fsrTime,mocapT,mocapR,mocapL,visualize)
% finds the locations from mocap 
% finds where the y-loc dips
% then confirms based on impacts found through fsr
% saves locations to coordinates matrix
% whichfoot: right is 1, left is 0
% 9/9/21, modified from before when it just used fsr to find the next coord
% now it inverts the y-sig, find peak, plus buffer, then store coord if ok
% with fsr
% extracted_idx_R/L stores the index of the extracted coordinate

impactT = fsrTime(impacts(:,1));
impactRorL = impacts(:,7);
coordinatesR = zeros(length(impactT),2);
coordinatesL = zeros(length(impactT),2);
coordinates = zeros(length(impactT),2);
whichfoot = zeros(length(impactT),1);
extracted_idx_R = zeros(length(impactT),1);
extracted_idx_L = zeros(length(impactT),1);

% right leg
mocapRavg = movmean(-1.*mocapR(:,2),10); % flip height sig to find "peaks"
[~,locsR,~,~] = findpeaks(mocapRavg,'MinPeakDistance',160); % messy data, so peaks are >70 pts apart
leg_still = 20; % leg stands still for at least 20 pts
locsR = locsR + leg_still;
% make sure it doesn't exceed mocapT length
locsRidx = locsR < length(mocapT);
locsRtime = mocapT(locsR(locsRidx));

% left leg
mocapLavg = movmean(-1.*mocapL(:,2),10); % flip height sig to find "peaks"
[~,locsL,~,~] = findpeaks(mocapLavg,'MinPeakDistance',160); % messy data, so peaks are >70 pts apart
locsL = locsL + leg_still;
% make sure it doesn't exceed mocapT length
locsLidx = locsL < length(mocapT);
locsLtime = mocapT(locsL(locsLidx));

for i = 1:length(impactT)
    heeltime = impactT(i);
    RorL = impactRorL(i);
    if RorL == 0 % left
        [~,index] = min(abs(locsLtime - heeltime));
        mocapI = locsL(index);
        coordinatesL(i,:) = [mocapL(mocapI,1),mocapL(mocapI,3)];
        coordinates(i,:) = [mocapL(mocapI,1),mocapL(mocapI,3)];
        whichfoot(i) = 0;
        extracted_idx_L(i) = mocapI;
    else
        [~,index] = min(abs(locsRtime - heeltime));
        mocapI = locsR(index);
        coordinatesR(i,1) = mocapR(mocapI,1);
        coordinatesR(i,2) = mocapR(mocapI,3);
        coordinates(i,:) = [mocapR(mocapI,1),mocapR(mocapI,3)];
        whichfoot(i) = 1;
        extracted_idx_R(i) = mocapI;
    end
end

% deleting 0 elements
Rzeros = find(extracted_idx_R == 0);
extracted_idx_R(Rzeros) = [];
coordinatesR(Rzeros,:) = [];
Lzeros = find(extracted_idx_L == 0);
extracted_idx_L(Lzeros) = [];
coordinatesL(Lzeros,:) = [];


% extracted_idx_R = nonzeros(extracted_idx_R);
% extracted_idx_L = nonzeros(extracted_idx_L);
% coordinatesR(~any(coordinatesR,2),:) = [];
% coordinatesL(~any(coordinatesL,2),:) = [];

% rightft = 1;
% avg_seg = 10; % take the average of 10 datapoints
% for i = 1:length(impactT)
%     starti = findTindex(impactT(i),mocapT);
%     endi = findTindex(impactT(i)+footstep_time,mocapT);
%     if isempty(endi)
%         endi = length(mocapT);
%     end
%     segL = mocapL(starti:endi,:);
%     segR = mocapR(starti:endi,:);
%     if std(segL(:,1)) > std(segR(:,1)) % left seg has more variation
%         % right leg is still
%         prev_avg = mean(segR(1:1+avg_seg,1)); % loop until x-value is steady
%         for j = avg_seg:(length(segR)-avg_seg)
%             curr_avg = mean(segR(j:j+avg_seg,1));
%             if abs(prev_avg-curr_avg) < 50 % movement is < 50 mm
%                 coordinatesR(i,1) = segR(j,1);
%                 coordinatesR(i,2) = segR(j,3); % store x,z coordinates
%                 whichfoot(i,1) = rightft;
%                 extracted_pts_R(i,1) = starti + j;
%                 coordinates(i,:) = [segR(j,1),segR(j,3)];
%                 break;
%             end
%             j = j+avg_seg;
%         end
%     else
%         % left leg is still
%         prev_avg = mean(segL(1:1+avg_seg,1)); % loop until x-value is steady
%         for j = avg_seg:(length(segL)-avg_seg)
%             curr_avg = mean(segL(j:j+avg_seg,1));
%             if abs(prev_avg-curr_avg) < 50 % movement is < 50 mm
%                 coordinatesL(i,1) = segL(j,1);
%                 coordinatesL(i,2) = segL(j,3); % store x,z coordinates
%                 coordinates(i,:) = [segL(j,1),segL(j,3)];
%                 extracted_pts_L(i,1) = starti + j;
%                 break;
%             end
%             j = j+avg_seg;
%         end
%     end
% end


if visualize
    figure;
    subplot(3,1,1)
    plot(mocapT,mocapR(:,1))
    hold on
    plot(mocapT(extracted_idx_R), coordinatesR(:,1),'r.','MarkerSize',10)
    title('Right foot x')
    subplot(3,1,2)
    plot(mocapT,mocapR(:,3))
    hold on
    plot(mocapT(extracted_idx_R), coordinatesR(:,2),'r.','MarkerSize',10)
    title('Right foot z')
    subplot(3,1,3)
    plot(mocapT,mocapR(:,2))
    hold on
    plot(mocapT(extracted_idx_R), mocapR(extracted_idx_R,2),'r.','MarkerSize',10)
    title('Right foot y')

    figure;
    subplot(3,1,1)
    plot(mocapT,mocapL(:,1))
    hold on
    plot(mocapT(extracted_idx_L), coordinatesL(:,1),'r.','MarkerSize',10)
    title('Left foot x')
    subplot(3,1,2)
    plot(mocapT,mocapL(:,3))
    hold on
    plot(mocapT(extracted_idx_L), coordinatesL(:,2),'r.','MarkerSize',10)
    title('Left foot z')
    subplot(3,1,3)
    plot(mocapT,mocapL(:,2))
    hold on
    plot(mocapT(extracted_idx_L), mocapL(extracted_idx_L,2),'r.','MarkerSize',10)
    title('Left foot y')
end


end

