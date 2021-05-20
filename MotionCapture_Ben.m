%% Motion capture data
motion_data = readtable('C:\Users\bpdwy\Dropbox (MIT)\Analysis\MotionCapture\04-20-2021-03-53-51.csv');

sensor = 12800;
motion = 180;

motion_first = 465;
sensor_first = 106750;

lx = table2array(motion_data(:,9));
ly = table2array(motion_data(:,11));
lz = table2array(motion_data(:,10));
lt = table2array(motion_data(:,2));

rx = table2array(motion_data(:,18));
ry = table2array(motion_data(:,20));
rz = table2array(motion_data(:,19));
rt = table2array(motion_data(:,2));

x = vertcat(lx, rx);
y = vertcat(ly, ry);
z = vertcat(lz, rz);
t = vertcat(lt, rt);

start_mag = 0;
start_idx = 0;

for i=1:length(x)
    if isnan(x(i))
        x(i) = 0;
    end
    if isnan(y(i))
        y(i) = 0;
    end
    if isnan(z(i))
        z(i) = 0;
    end
    if z(i) > start_mag
        start_mag = z(i);
        start_idx = i - length(x)/2;
    end
    if z(i) > 350
        z(i) = 350;
    end
    if isnan(t(i))
        t(i) = 0;
    end
end

order = [];
[peaks_l, idx_l] = findpeaks(-lz(1:2200),'MinPeakDistance', 200);
[peaks_r, idx_r] = findpeaks(-rz(1:2200),'MinPeakDistance', 200);
step_l = zeros(length(idx_l),3);
step_r = zeros(length(idx_r),3);
for s = 1:length(idx_l)
    step_l_idx = idx_l(s);
    step_r_idx = idx_r(s);
    step_l(s, :) = [lx(step_l_idx) -ly(step_l_idx) step_l_idx];
    step_r(s, :) = [rx(step_l_idx) -ry(step_l_idx) step_r_idx];
end

stomp_idx = 0;
for num = length(step_r):-1:1
    if step_r(num,3) > start_idx
        stomp_idx = step_r(num, 3);
    end
end

sensor_x = [-3585 -3606 3616 3607];
sensor_y = [3395 -196 -212 3376];

joined_steps = zeros(length(step_r)*2,3);
spot = 1;
for i = 1:length(step_r)
    joined_steps(spot, :) = step_r(i,:);
    spot = spot + 1;
    joined_steps(spot, :) = step_l(i,:);
    spot = spot + 1;
end

dst_1 = zeros(length(joined_steps),1);
dst_2 = zeros(length(joined_steps),1);
dst_3 = zeros(length(joined_steps),1);
dst_4 = zeros(length(joined_steps),1);

for r = 1:length(joined_steps)
    dst_1(r) = power((power((sensor_x(1) - joined_steps(r, 1)),2) + power((sensor_y(1) - joined_steps(r, 2)),2)),1/2);
    dst_2(r) = power((power((sensor_x(2) - joined_steps(r, 1)),2) + power((sensor_y(2) - joined_steps(r, 2)),2)),1/2);
    dst_3(r) = power((power((sensor_x(3) - joined_steps(r, 1)),2) + power((sensor_y(3) - joined_steps(r, 2)),2)),1/2);
    dst_4(r) = power((power((sensor_x(4) - joined_steps(r, 1)),2) + power((sensor_y(4) - joined_steps(r, 2)),2)),1/2);
end

dst = horzcat(dst_1, dst_2, dst_3, dst_4);

%% Plots
figure();
hold on
scatter(x(1:2200),-y(1:2200),10, z(1:2200));
scatter(step_l(:,1), step_l(:, 2), 75, 'MarkerEdgeColor', 'magenta',... 
                                        'MarkerFaceColor', 'none',...
                                        'LineWidth', 2);
xlabel('x');
ylabel('y');
c = colorbar;
c.Label.String = 'z';
grid on;

Right 
figure();
hold on
scatter(x(3254:5454),-y(3254:5454),10, z(3254:5454));
scatter(step_r(:,1), step_r(:, 2), 75, 'MarkerEdgeColor', 'red',... 
                                        'MarkerFaceColor', 'none',...
                                        'LineWidth', 2);
scatter(sensor_x, sensor_y, 75, 'MarkerEdgeColor', 'black',... 
                                'MarkerFaceColor', 'black',...
                                'LineWidth', 2);
xlabel('x');
ylabel('y');
c = colorbar;
c.Label.String = 'z';
grid on;

%Both
figure();
hold on
left = scatter(x(1:2200),-y(1:2200),10, z(1:2200));
right = scatter(x(3254:5454),-y(3254:5454),10, z(3254:5454));
scatter(step_l(:,1), step_l(:, 2), 75, 'MarkerEdgeColor', 'magenta',... 
                                        'MarkerFaceColor', 'none',...
                                        'LineWidth', 2);
scatter(step_r(:,1), step_r(:, 2), 75, 'MarkerEdgeColor', 'red',... 
                                        'MarkerFaceColor', 'none',...
                                        'LineWidth', 2);
scatter(sensor_x, sensor_y, 75, 'MarkerEdgeColor', 'black',... 
                                'MarkerFaceColor', 'black',...
                                'LineWidth', 2);
xlabel('x');
ylabel('y');
c = colorbar;
c.Label.String = 'z';
grid on;


%% Sensor 
sensor_data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\ImmersionMatlab\data_04-20-2021_15-54-17.mat');
sensor_sum = sum(abs(sensor_data.datas), 2);

[peaks1, idx1, width1, prominence1] = findpeaks(sensor_sum(100000:352000),... 
                                                'MinPeakHeight', 0.025,... 
                                                'MinPeakDistance', 5000,... 
                                                'MinPeakProminence', 0.025);
                                            
time = sensor_data.times;
filt_time = time - time(100000);
final_time = filt_time(100000:352000);

peak_times = idx1/12800;

%% S
sensor_data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\ImmersionFloorSensors\data_04-20-2021_15-54-17.mat').datas;
max_val = 0;
max_idx = 0;
for i = 1:length(sensor_data)
    for c = 1:4
        if sensor_data(i,c) > max_val
            max_val = sensor_data(i,c);
            max_idx = i;
        end
    end
end

%% Sensor distances
sensor_data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\ImmersionFloorSensors\data_04-20-2021_15-54-17.mat').datas;

[peaks1, idx1, width1, prominence1] = findpeaks(sensor_data(100000:352000,1),... 
                                                'MinPeakHeight', 0.005,... 
                                                'MinPeakDistance', 2500,... 
                                                'MinPeakProminence', 0.005);
[peaks2, idx2, width2, prominence2] = findpeaks(sensor_data(100000:352000,2),... 
                                                'MinPeakHeight', 0.005,... 
                                                'MinPeakDistance', 2500,... 
                                                'MinPeakProminence', 0.005);
[peaks3, idx3, width3, prominence3] = findpeaks(sensor_data(100000:352000,3),... 
                                                'MinPeakHeight', 0.005,... 
                                                'MinPeakDistance', 2500,... 
                                                'MinPeakProminence', 0.005);
[peaks4, idx4, width4, prominence4] = findpeaks(sensor_data(100000:352000,4),... 
                                                'MinPeakHeight', 0.005,... 
                                                'MinPeakDistance', 2500,... 
                                                'MinPeakProminence', 0.005);
                                            
figure();
hold on
% plot(peaks1);
% plot(peaks2);
plot(peaks3);
plot(peaks4);

figure();
hold on
% plot(dst(:, 1));
% plot(dst(:, 2));
plot(dst(:, 3));
plot(dst(:, 4));
%% Time Shift & Peaks
sensor_data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\ImmersionFloorSensors\data_04-20-2021_15-54-17.mat').datas;
motion_data = readtable('C:\Users\bpdwy\Dropbox (MIT)\Analysis\MotionCapture\04-20-2021-03-53-51.csv');

[motion_data_crop,sensor_data_crop, sensor_start, motion_start] = TimeShift(motion_data,sensor_data);

[peaks1, peaks2, peaks3, peaks4, idx1, idx2, idx3, idx4] = TDOA4new(sensor_data_crop);

%% MoCap func
sensor_data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\ImmersionFloorSensors\data_04-20-2021_15-54-17.mat').datas;
motion_data = readtable('C:\Users\bpdwy\Dropbox (MIT)\Analysis\MotionCapture\04-20-2021-03-53-51.csv');

[motion_data_crop,sensor_data_crop, sensor_start, motion_start] = TimeShift(motion_data,sensor_data);

[joined_steps, dst] = ProcessMoCap(motion_data_crop);
[sensor_data_crop_clean] = lpf_data(sensor_data_crop);
[onset_idx, peak_idx, peak_val] = TDOA_edited(sensor_data_crop,length(joined_steps));

adjusted = zeros(size(peak_idx));
sen_1 = horzcat(peak_idx(1,:).', peak_val(1,:).');
sen_2 = horzcat(peak_idx(2,:).', peak_val(2,:).');
sen_3 = horzcat(peak_idx(3,:).', peak_val(3,:).');
sen_4 = horzcat(peak_idx(4,:).', peak_val(4,:).');

sen_1_sorted = sortrows(sen_1, 1);
sen_2_sorted = sortrows(sen_2, 1);
sen_3_sorted = sortrows(sen_3, 1);
sen_4_sorted = sortrows(sen_4, 1);

sen_merged_1 = horzcat(sen_1_sorted(:, 1), sen_2_sorted(:, 1), sen_3_sorted(:, 1), sen_4_sorted(:, 1));
sen_merged_2 = horzcat(sen_1_sorted(:, 2), sen_2_sorted(:, 2), sen_3_sorted(:, 2), sen_4_sorted(:, 2));
sen_merged = horzcat(sen_merged_1, sen_merged_2);

final = zeros(length(joined_steps), 12);
for r = 1:length(joined_steps)
    first_idx = 1000000000000;
    high_mag = 0;
    for c = 1:4
        if sen_merged(r, c) < first_idx
            first_idx = sen_merged(r, c);
        end
    end
    for c = 5:8
        if sen_merged(r, c) > high_mag
            high_mag = sen_merged(r, c);
        end
    end
    for c = 1:4
        final(r, c) = sen_merged(r, c) - first_idx;
    end
    for c = 5:8
        final(r, c) = sen_merged(r, c) - high_mag;
    end
    for c = 9:12
        final(r,c) = dst(r,c-8);
    end
end

%% Combining MoCap

sensor_data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\ImmersionFloorSensors\data_04-20-2021_15-54-17.mat').datas;
motion_data = readtable('C:\Users\bpdwy\Dropbox (MIT)\Analysis\MotionCapture\04-20-2021-03-53-51.csv');
[final] = Prep_MoCap(motion_data, sensor_data);