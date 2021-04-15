%% 3/12/2021
depressed = readtable("C:/Users/bpdwy/Documents/MIT/Academic/UROP/Depressed.csv");
control = readtable("C:/Users/bpdwy/Documents/MIT/Academic/UROP/Control.csv");
combo = vertcat(depressed(:,2:7), control(:,2:7));
labels = combo(:,2);
normalized = zeros(size(combo, 1),size(combo, 2));
for i = 1:6
    for h = 1:200
        if i == 2
            normalized(h,i) = combo{h,i};
        end
        minimum = min(combo{:,i});
        maximum = max(combo{:,i});
        normalized(h, i) = (combo{h,i} - minimum)/(maximum - minimum);
    end
end

%% 10 percent reduced cadence
depressed = readtable("C:/Users/bpdwy/Documents/MIT/Academic/UROP/Depressed_90_cadence.csv");
control = readtable("C:/Users/bpdwy/Documents/MIT/Academic/UROP/Control_90_cadence.csv");
combo = vertcat(depressed(:,2:7), control(:,2:7));
labels = combo(:,2);
normalized = zeros(size(combo, 1),size(combo, 2));
for i = 1:6
    for h = 1:200
        if i == 2
            normalized(h,i) = combo{h,i};
        end
        minimum = min(combo{:,i});
        maximum = max(combo{:,i});
        normalized(h, i) = (combo{h,i} - minimum)/(maximum - minimum);
    end
end

%% 25 percent reduced cadence
depressed = readtable("C:/Users/bpdwy/Documents/MIT/Academic/UROP/Depressed_75_cadence.csv");
control = readtable("C:/Users/bpdwy/Documents/MIT/Academic/UROP/Control_75_cadence.csv");
combo = vertcat(depressed(:,2:7), control(:,2:7));
labels = combo(:,2);
normalized = zeros(size(combo, 1),size(combo, 2));
for i = 1:6
    for h = 1:200
        if i == 2
            normalized(h,i) = combo{h,i};
        end
        minimum = min(combo{:,i});
        maximum = max(combo{:,i});
        normalized(h, i) = (combo{h,i} - minimum)/(maximum - minimum);
    end
end

%% Noise 0-1
depressed = readtable("C:/Users/bpdwy/Documents/MIT/Academic/UROP/Depressed_noise_0_1.csv");
control = readtable("C:/Users/bpdwy/Documents/MIT/Academic/UROP/Control_noise_0_1.csv");
combo = vertcat(depressed(:,2:7), control(:,2:7));
combo.StepWidth = [];
combo.Stride = combo.Stride/2;
labels = combo(:,2);
normalized = zeros(size(combo, 1),size(combo, 2));
for i = 1:5
    for h = 1:200
        if i == 2
            normalized(h,i) = combo{h,i};
        end
        minimum = min(combo{:,i});
        maximum = max(combo{:,i});
        normalized(h, i) = (combo{h,i} - minimum)/(maximum - minimum);
    end
end

%% Noise 0-3
depressed = readtable("C:/Users/bpdwy/Documents/MIT/Academic/UROP/Depressed_noise_0_3.csv");
control = readtable("C:/Users/bpdwy/Documents/MIT/Academic/UROP/Control_noise_0_3.csv");
combo = vertcat(depressed(:,2:7), control(:,2:7));
combo.StepWidth = [];
combo.Stride = combo.Stride/2;
labels = combo(:,2);
normalized = zeros(size(combo, 1),size(combo, 2));
for i = 1:5
    for h = 1:200
        if i == 2
            normalized(h,i) = combo{h,i};
        end
        minimum = min(combo{:,i});
        maximum = max(combo{:,i});
        normalized(h, i) = (combo{h,i} - minimum)/(maximum - minimum);
    end
end

%% Noise mean
depressed = readtable("C:/Users/bpdwy/Documents/MIT/Academic/UROP/Depressed_noise_mean_17.csv");
control = readtable("C:/Users/bpdwy/Documents/MIT/Academic/UROP/Control_noise_mean_17.csv");
combo = vertcat(depressed(:,2:7), control(:,2:7));
combo.StepWidth = [];
combo.Stride = combo.Stride/2;
labels = combo(:,2);
normalized = zeros(size(combo, 1),size(combo, 2));
for i = 1:5
    for h = 1:200
        if i == 2
            normalized(h,i) = combo{h,i};
        end
        minimum = min(combo{:,i});
        maximum = max(combo{:,i});
        normalized(h, i) = (combo{h,i} - minimum)/(maximum - minimum);
    end
end

%% Mean and SD of our data

steps = load('C:\Users\bpdwy\Documents\MIT\Academic\UROP\lab_analysis\data\data_10-19-2020_15-56.mat').datas;
[peaks1, peaks2, peaks3, idx1, idx2, idx3] = TDOA4new(steps);
num_steps = length(peaks1);
time_between = [];
Fs = 12800;
for i = 1:num_steps-1
%     first_first = min([idx1(i), idx2(i), idx3(i)]);
%     second_first = min([idx1(i+1), idx2(i+1), idx3(i+1)]);
    val = (idx1(i+1) - idx1(i))/Fs;
    time_between = [time_between, val]; 
end
avg_time = mean(time_between); %sec/step
step_per_sec = 1/avg_time;
cadence = 60*step_per_sec;

%% test function indoor walking
data = struct("data1",load('C:\Users\bpdwy\Documents\MIT\Academic\UROP\lab_analysis\data\data_10-19-2020_15-49.mat').datas, ...
"data2",load('C:\Users\bpdwy\Documents\MIT\Academic\UROP\lab_analysis\data\data_10-19-2020_15-56.mat').datas,...
"data3",load('C:\Users\bpdwy\Documents\MIT\Academic\UROP\lab_analysis\data\data_10-19-2020_15-58.mat').datas,...
"data4",load('C:\Users\bpdwy\Documents\MIT\Academic\UROP\lab_analysis\data\data_10-19-2020_15-59.mat').datas,...
"data5",load('C:\Users\bpdwy\Documents\MIT\Academic\UROP\lab_analysis\data\data_10-19-2020_16-00.mat').datas,...
"data6",load('C:\Users\bpdwy\Documents\MIT\Academic\UROP\lab_analysis\data\data_10-19-2020_16-01.mat').datas,...
"data7",load('C:\Users\bpdwy\Documents\MIT\Academic\UROP\lab_analysis\data\data_10-19-2020_16-02.mat').datas,...
"data8",load('C:\Users\bpdwy\Documents\MIT\Academic\UROP\lab_analysis\data\data_10-19-2020_16-03.mat').datas,...
"data9",load('C:\Users\bpdwy\Documents\MIT\Academic\UROP\lab_analysis\data\data_10-19-2020_16-04.mat').datas,...
"data10",load('C:\Users\bpdwy\Documents\MIT\Academic\UROP\lab_analysis\data\data_10-19-2020_16-05.mat').datas);
[cadence, stan, vel] = WalkingCadence(data);
