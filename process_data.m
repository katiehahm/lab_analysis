% 1/24/21
% used to run save_processed_data fn to record features into loadable workspace file
% ran on data_12_10_2020 to 01_15_2021


filename = 'data_01-15-2021_12-13-40';
n = 10;
sn = 4;
m = 0.26;
save_processed_data(filename, n, sn, m)

%%
% 1/24/21
% to delete noise that interfers with peak finding alg

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\data\';
load([data_root_katie, filename])
minidx = 795700;
maxidx = 839800;
datas(minidx:maxidx,:) = 0;
times(minidx:maxidx,:) = 0;
save([data_root_katie, filename], 'datas', 'times', '-append')