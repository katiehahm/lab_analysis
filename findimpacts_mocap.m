function [extracted_pts_R, extracted_pts_L, coordinates, whichfoot] = findimpacts_mocap(impactT,mocapT,mocapR,mocapL,visualize)
% finds the locations from mocap based on impacts found through fsr
% saves locations to coordinates matrix
% 6/8/21

footstep_time = 0.3; % length of time foot is reliably still
coordinatesR = ones(length(impactT),2);
coordinatesL = ones(length(impactT),2);
coordinates = zeros(length(impactT),2);
whichfoot = zeros(length(impactT),1);
extracted_pts_R = ones(length(impactT),1);
extracted_pts_L = ones(length(impactT),1);
rightft = 1;
leftft = 2;
avg_seg = 10; % take the average of 10 datapoints
for i = 1:length(impactT)
    starti = findTindex(impactT(i),mocapT);
    endi = findTindex(impactT(i)+footstep_time,mocapT);
    if isempty(endi)
        endi = length(mocapT);
    end
    segL = mocapL(starti:endi,:);
    segR = mocapR(starti:endi,:);
    if std(segL(:,1)) > std(segR(:,1)) % left seg has more variation
        % right leg is still
        prev_avg = mean(segR(1:1+avg_seg,1)); % loop until x-value is steady
        for j = avg_seg:(length(segR)-avg_seg)
            curr_avg = mean(segR(j:j+avg_seg,1));
            if abs(prev_avg-curr_avg) < 50 % movement is < 50 mm
                coordinatesR(i,1) = segR(j,1);
                coordinatesR(i,2) = segR(j,3); % store x,z coordinates
                whichfoot(i,1) = rightft;
                extracted_pts_R(i,1) = starti + j;
                coordinates(i,:) = [segR(j,1),segR(j,3)];
                break;
            end
            j = j+avg_seg;
        end
    else
        % left leg is still
        prev_avg = mean(segL(1:1+avg_seg,1)); % loop until x-value is steady
        for j = avg_seg:(length(segL)-avg_seg)
            curr_avg = mean(segL(j:j+avg_seg,1));
            if abs(prev_avg-curr_avg) < 50 % movement is < 50 mm
                coordinatesL(i,1) = segL(j,1);
                coordinatesL(i,2) = segL(j,3); % store x,z coordinates
                coordinates(i,:) = [segL(j,1),segL(j,3)];
                extracted_pts_L(i,1) = starti + j;
                break;
            end
            j = j+avg_seg;
        end
    end
end


if visualize
    figure;
    subplot(3,1,1)
    plot(mocapT,mocapR(:,1))
    hold on
    plot(mocapT(extracted_pts_R), coordinatesR(:,1),'r.','MarkerSize',10)
    title('Right foot x')
    subplot(3,1,2)
    plot(mocapT,mocapR(:,3))
    hold on
    plot(mocapT(extracted_pts_R), coordinatesR(:,2),'r.','MarkerSize',10)
    title('Right foot z')
    subplot(3,1,3)
    plot(mocapT,mocapR(:,2))
    hold on
    plot(mocapT(extracted_pts_R), mocapR(extracted_pts_R,2),'r.','MarkerSize',10)
    title('Right foot y')

    figure;
    subplot(3,1,1)
    plot(mocapT,mocapL(:,1))
    hold on
    plot(mocapT(extracted_pts_L), coordinatesL(:,1),'r.','MarkerSize',10)
    title('Left foot x')
    subplot(3,1,2)
    plot(mocapT,mocapL(:,3))
    hold on
    plot(mocapT(extracted_pts_L), coordinatesL(:,2),'r.','MarkerSize',10)
    title('Left foot z')
    subplot(3,1,3)
    plot(mocapT,mocapL(:,2))
    hold on
    plot(mocapT(extracted_pts_L), mocapL(extracted_pts_L,2),'r.','MarkerSize',10)
    title('Left foot y')
end

size(coordinates)

end
