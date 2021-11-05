% 9/9/21
% uses left and right directions of mocap to extract out turning points
% look at plots of distance changes to set a threshold for where turns are

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj';
subj = '2'; % number of subject
intervention = 'normal1'; % normal1 slow insole weight count normal2
filename = [data_root_katie, subj, '_', intervention];

load(filename)

%% 10/20/21 
% assumes that subj always walks in the long direction, so 
% any change in x direction is a turn

% compare these two figures with the two at the end of this section
% the turning points should've been eliminated
figure;
plot(arrival_idx(:,1),coordinates(:,1),'.')
title('All footfalls')
xlabel('Sample number')
ylabel('x coordinates')
figure;
plot(arrival_idx(:,1),coordinates(:,2),'.')


data = zeros(length(impacts),19);
new_impacts = zeros(length(impacts),7);
% start of segment is -1, end of segment is 1, 0 in between
toe_t = fsrTime(impacts(1,4));
heel_t = fsrTime(impacts(1,1));
data(1,:) = [arrival_idx(1,:),peak_idx(1,:),peak_mag(1,:),coordinates(1,:),...
    impacts(1,3),toe_t-heel_t, NaN, NaN, -1];
new_impacts(1,:) = impacts(1,:);
curr_x = sign( coordinates(2,1)-coordinates(1,1) );
index_count = 2;
for i = 2:(length(coordinates)-1)
    new_x = sign( coordinates(i+1,1)-coordinates(i,1) );
    if curr_x ~= new_x % change in direction
        curr_x = new_x;
        % delete datapoint and end/start the segment
        data(index_count-1,19) = 1;
        data(index_count,19) = -1;
    else % no change in direction
        toe_t = fsrTime(impacts(i+1,4));
        heel_t = fsrTime(impacts(i+1,1));
        data(index_count,1:18) = [arrival_idx(i+1,:),peak_idx(i+1,:),peak_mag(i+1,:),coordinates(i+1,:),...
            impacts(i+1,3),toe_t-heel_t, coordinates(i,:)];
        index_count = index_count+1;
        new_impacts(index_count,:) = impacts(i+1,:);
    end
end

% delete the empty rows at the end of data
data(~any(data,2),:) = [];

figure;
plot(data(:,1),data(:,13),'.')
title('Only straight path footfalls')
xlabel('Sample number')
ylabel('x coordinates')
figure;
plot(data(:,1),data(:,14),'.')

% now run the save file section at the end

%% 9/9/21
% this and next section relies on step length inconsistencies
% to parse out where the turns are
coordN = length(whichfoot);
coordinatesR = zeros(length(extracted_pts_R),2);
coordinatesL = zeros(length(extracted_pts_L),2);
Rcount = 1;
Lcount = 1;

for i = 1:length(whichfoot)
    if whichfoot(i) == 1
        coordinatesR(Rcount,:) = coordinates(i,:);
        Rcount = Rcount + 1;
    else
        coordinatesL(Lcount,:) = coordinates(i,:);
        Lcount = Lcount + 1;
    end
end

% right foot
% figure;
% plot(mocapT,mocapR(:,1))
% hold on
% plot(mocapT(extracted_pts_R),mocapR(extracted_pts_R,1),'ro')
% 
% step_lengthR = zeros(1,length(coordinatesR)-1);
% for i = 2:length(coordinatesR)
%     distance = distance_btw_coordinates(extracted_pts_R(i-1),extracted_pts_R(i),mocapR, 1);
%     step_lengthR(i-1) = distance;
% end
% 
% figure; plot(mocapT(extracted_pts_R(2:end,:)),step_lengthR, '.')
% ylim([0 3])
% 
% % left foot
% figure;
% plot(mocapT,mocapL(:,1))
% hold on
% plot(mocapT(extracted_pts_L),mocapL(extracted_pts_L,1),'ro')
% 
% step_lengthL = zeros(1,length(coordinatesL)-1);
% for i = 2:length(coordinatesL)
%     distance = distance_btw_coordinates(extracted_pts_L(i-1),extracted_pts_L(i),mocapL, 1);
%     step_lengthL(i-1) = distance;
% end
% 
% figure; plot(mocapT(extracted_pts_L(2:end,:)),step_lengthL, '.')
% ylim([0 3])
% 
% %% extract turns (x-dir)
% 
% % set these values based on plots above
% thresh_min = 0.895;
% thresh_max = 1.8;


figure; 
plot(mocapT(extracted_pts_L),mocapL(extracted_pts_L,1),'r.')
hold on
plot(mocapT(extracted_pts_R),mocapR(extracted_pts_R,1),'b.')
step_lengths = zeros(1,length(coordinates)-1);
coordinates_idx = [extracted_pts_L; extracted_pts_R];
coordinates_idx = sort(coordinates_idx);
for i = 2:length(coordinates)
    x_curr = coordinates(i,1);
    y_curr = coordinates(i,2);
    x_prev = coordinates(i-1,1);
    y_prev = coordinates(i-1,2);
    if isnan(x_curr) || isnan(x_prev)
        distance = abs(y_curr-y_prev);
    elseif isnan(y_curr) || isnan(y_prev)
        distance = abs(x_curr-x_prev);
    else
        distance = sqrt( (x_curr-x_prev)^2 + (y_curr-y_prev)^2 );
    end
    step_lengths(i-1) = distance;
end
figure;
plot(mocapT(coordinates_idx(2:end)),step_lengths,'k.')

%% 9/9/21 continuation: extract turns 

% set these values based on plots above
thresh_min = 0.379;
thresh_max = 1;

data = zeros(length(impacts),19);
% start of segment is -1, end of segment is 1, 0 in between
toe_t = fsrTime(impacts(1,4));
heel_t = fsrTime(impacts(1,1));
data(1,:) = [arrival_idx(1,:),peak_idx(1,:),peak_mag(1,:),coordinates(1,:),...
    impacts(1,3),toe_t-heel_t, NaN, NaN, -1];
index_count = 2;
for i = 1:length(step_lengths)
    if step_lengths(i) > thresh_min && step_lengths(i) < thresh_max % ok
        toe_t = fsrTime(impacts(i+1,4));
        heel_t = fsrTime(impacts(i+1,1));
        data(index_count,1:18) = [arrival_idx(i+1,:),peak_idx(i+1,:),peak_mag(i+1,:),coordinates(i+1,:),...
            impacts(i+1,3),toe_t-heel_t, coordinates(i,:)];
        index_count = index_count+1;
    else
        % delete datapoint and end/start the segment
        data(index_count-1,19) = 1;
        data(index_count,19) = -1;
    end
end

% delete the empty rows at the end of data
data(~any(data,2),:) = [];

%% save data matrix to excel and matlab workspace
filename = [filename, '_extract_straight_paths'];
% save to matlab
save(filename,'data','new_impacts')

% save to excel
excel_filename = [filename,'_excel.xlsx'];
T = table(data);
writetable(T,excel_filename);

