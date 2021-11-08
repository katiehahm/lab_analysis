function [extracted_pts_R,extracted_pts_L,coordinates] = manual_fix_mocap(r_wrong,r_right,l_wrong,l_right,mocapR,mocapL,extracted_pts_R,extracted_pts_L,coordinates,whichfoot)
% takes in manual numbers to change for mocap coordinates extraction
% inputs are arrays
% 11/5/21

for i = 1:length(r_wrong)
    indeces = abs(extracted_pts_R - r_wrong(i));
    [~,change_idx] = min(indeces);
    extracted_pts_R(change_idx) = r_right(i);
    
    coord_count = 0;
    coord_change_idx = 0;
    for j = 1:length(coordinates)
        if whichfoot(j) == 1 % right foot
            coord_count = coord_count + 1;
            if coord_count == change_idx
                coord_change_idx = j;
                break;
            end
        end
    end
    coordinates(coord_change_idx,:) = [mocapR(extracted_pts_R(change_idx),1),mocapR(extracted_pts_R(change_idx),3)];
end

for i = 1:length(l_wrong)
    indeces = abs(extracted_pts_L - l_wrong(i));
    [~,change_idx] = min(indeces);
    extracted_pts_L(change_idx) = l_right(i);
    
    coord_count = 0;
    coord_change_idx = 0;
    for j = 1:length(coordinates)
        if whichfoot(j) == 0 % left foot
            coord_count = coord_count + 1;
            if coord_count == change_idx
                coord_change_idx = j;
                break;
            end
        end
    end
    coordinates(coord_change_idx,:) = [mocapL(extracted_pts_L(change_idx),1),mocapL(extracted_pts_L(change_idx),3)];
end

% visually check
figure;
subplot(3,1,1)
plot(mocapR(:,1))
hold on
plot(extracted_pts_R, mocapR(extracted_pts_R,1),'r.','MarkerSize',10)
title('Right foot x')
subplot(3,1,2)
plot(mocapR(:,3))
hold on
plot(extracted_pts_R, mocapR(extracted_pts_R,3),'r.','MarkerSize',10)
title('Right foot z')
subplot(3,1,3)
plot(mocapR(:,2))
hold on
plot(extracted_pts_R, mocapR(extracted_pts_R,2),'r.','MarkerSize',10)
title('Right foot y')

figure;
subplot(3,1,1)
plot(mocapL(:,1))
hold on
plot(extracted_pts_L, mocapL(extracted_pts_L,1),'r.','MarkerSize',10)
title('Left foot x')
subplot(3,1,2)
plot(mocapL(:,3))
hold on
plot(extracted_pts_L, mocapL(extracted_pts_L,3),'r.','MarkerSize',10)
title('Left foot z')
subplot(3,1,3)
plot(mocapL(:,2))
hold on
plot(extracted_pts_L, mocapL(extracted_pts_L,2),'r.','MarkerSize',10)
title('Left foot y')
end