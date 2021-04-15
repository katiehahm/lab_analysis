%% Dummy Data
throw_1_data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_01-25-2021_17-25-02').datas;
throw_2_data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_01-25-2021_17-25-22').datas;
throw_3_data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_01-25-2021_17-25-58').datas;

drop_1_data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_01-25-2021_15-27-13').datas;
drop_2_data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_01-25-2021_15-34-42').datas;

steps_1_data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_10-02-2020_19-03').datas;

%% Creating an "impact_length" feature
throw_1_data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_01-25-2021_17-25-02').datas;
throw_2_data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_01-25-2021_17-25-22').datas;
throw_3_data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_01-25-2021_17-25-58').datas;

drop_1_data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_01-25-2021_15-27-13').datas;
drop_2_data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_01-25-2021_15-34-42').datas;

step_1_data = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_10-02-2020_19-03').datas(138000:144000, :);

sensor_sum = sum(abs(step_1_data), 2);
final = movmean(sensor_sum, 250);
start = "True";
for dat = 1:size(final, 1)
    if start == "True"
        if final(dat) > 0.05
            if start == "True"
                begin = dat;
                start = "False";
            end
        end
    else
        if final(dat) < 0.05
            stop = dat;
        end
    end
end
impact_length = stop - begin;

% drop 1 length = 1722
% drop 2 length = 1591
% throw 1 length = 2395
% throw 2 length = 2873
% throw 3 length = 1653
% step 1 length = 