%% 6/15/21 form one csv file for python code

clear all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\ProcessedData\';

datestr = '06-04-2021_13-47-54';
load([data_root_katie, datestr])
DataT(any(isnan(DataT),2),:) = [];
final_data = DataT;

datestr = '06-10-2021_13-23-05';
load([data_root_katie, datestr])
DataT(any(isnan(DataT),2),:) = [];
final_data = [final_data;DataT];

datestr = '06-10-2021_13-38-16';
load([data_root_katie, datestr])
DataT(any(isnan(DataT),2),:) = [];
final_data = [final_data;DataT];

datestr = '06-10-2021_14-00-15';
load([data_root_katie, datestr])
DataT(any(isnan(DataT),2),:) = [];
final_data = [final_data;DataT];

datestr = '06-10-2021_14-17-29';
load([data_root_katie, datestr])
DataT(any(isnan(DataT),2),:) = [];
final_data = [final_data;DataT];

writematrix(final_data,[data_root_katie,'06-04_06-10_finaldata_nonan.csv'])

%% 6/18/21 fixing data csv to put relative index values instead of absolute

clear all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\ProcessedData\';
T = readtable([data_root_katie,'06-04_06-10_finaldata_nonan.csv']);
data = table2array(T);
arrival_idx = data(:,1:4);
peak_idx = data(:,5:8);
nrows = length(arrival_idx);
new_arrival = zeros(nrows,4);
new_peak = zeros(nrows,4);
for i = 1:nrows
    firstidx = min(arrival_idx(i,:));
    new_arrival(i,:) = arrival_idx(i,:) - firstidx;
    new_peak(i,:) = peak_idx(i,:) - arrival_idx(i,:);
end

new_final_data = [new_arrival, new_peak, data(:,9:end)];
writematrix(new_final_data,[data_root_katie,'06-04_06-10_edited_finaldata_nonan.csv'])

%% 6/21/21 changing coordinate labels to grid labels 1-18 for unsupervised

clear all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\ProcessedData\';
T = readtable([data_root_katie,'06-04_06-10_edited_finaldata_nonan.csv']);
data = table2array(T);
coords = data(:,13:end);

% adjust so that S1 is at (0,0)
coords = coords + [3584, 3394];

zmax = 3394 + 196;
xmax = 3584 + 3607;

nrows = length(coords);
grid_label = zeros(nrows,1);

zmod = round(zmax/3);
xmod = round(xmax/6);

% label each coordinate 1-18, S1 is 1, S2 is 3, S4 is 16, S3 is 18
for i = 1:nrows
    coord_x = coords(i,1);
    coord_z = coords(i,2);
    znum = floor(coord_z/zmod);
    xnum = floor(coord_x/xmod);
    label = 3*xnum + znum + 1;
    if label < 1
        label = 1;
    elseif label > 18
        label = 18;
    end
    grid_label(i) = label;
end

new_data = [data(:,1:12),grid_label];
writematrix(new_data,[data_root_katie,'06-04_06-10_edited_finaldata_nonan_gridlabel.csv'])
