%% 31.5" Black
data1 = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_12-10-2020_12-16-21').peak_mag;
data2 = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_12-10-2020_12-18-11').peak_mag;
data3 = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_12-10-2020_12-19-54').peak_mag;
data4 = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_12-10-2020_12-21-17').peak_mag;

datas = {data1, data2, data3, data4};
codes = [11100000000 00011100000 00000011100 00000000011];
[sensor_1, sensor_2, sensor_3, sensor_4] = combine_multiple_paths_ben(datas, codes);

sens1_coeffs = polyfit(sensor_1(:,2), sensor_1(:,1), 2)
sens2_coeffs = polyfit(sensor_2(:,2), sensor_2(:,1), 2)
sens3_coeffs = polyfit(sensor_3(:,2), sensor_3(:,1), 2)
sens4_coeffs = polyfit(sensor_4(:,2), sensor_4(:,1), 2)

%% 17" Black
data1 = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_01-05-2021_13-28-34').peak_mag;
data2 = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_01-05-2021_13-29-59').peak_mag;
data3 = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_01-05-2021_13-31-27').peak_mag;
data4 = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_01-05-2021_13-32-32').peak_mag;

datas = {data1, data2, data3, data4};
codes = [11100000000 00011100000 00000011100 00000000011];
[sensor_1, sensor_2, sensor_3, sensor_4] = combine_multiple_paths_ben(datas, codes);

sens1_coeffs = polyfit(sensor_1(:,2), sensor_1(:,1), 2)
sens2_coeffs = polyfit(sensor_2(:,2), sensor_2(:,1), 2)
sens3_coeffs = polyfit(sensor_3(:,2), sensor_3(:,1), 2)
sens4_coeffs = polyfit(sensor_4(:,2), sensor_4(:,1), 2)

%% 31.5" White 
data1 = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_01-25-2021_15-27-13').peak_mag;
data2 = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_01-05-2021_13-38-22').peak_mag;
data3 = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_01-05-2021_13-40-13').peak_mag;
data4 = load('C:\Users\bpdwy\Dropbox (MIT)\Analysis\data\data_01-05-2021_13-41-42').peak_mag;

datas = {data1, data2, data3, data4};
codes = [11100000000 00011100000 00000011100 00000000011];
[sensor_1, sensor_2, sensor_3, sensor_4] = combine_multiple_paths_ben(datas, codes);