% to make excel for localization 11/30/21

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\';
takes = {'regular1', 'regular2', 'slow', 'stiffRknee', 'stiffLknee', 'insoleRweightL1', 'insoleRweightL2'};

featureV = zeros(1,22);
for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie, intervention];
    load(string(filename))
    
    for i = 2:length(impacts)
        if walk_segments(i) ~= -1 % don't count prev step if it's a new walking segment
            xcoord = coordinates(i,1)/1000;
            prev_x = coordinates(i-1,1)/1000;
            
            curr_mag = peak_mag(i,:);
            prev_mag = peak_mag(i-1,:);
            curr_energy = energy_squared(i,:);
            prev_energy = energy_squared(i-1,:);
            
            min_arr = min(arrival_idx(i,:));
            arrivals = arrival_idx(i,:) - min_arr;
            
            % xy coord, arrival idx, peak mag, energy, xy prev coord
            feature = [xcoord,curr_mag, arrivals, curr_energy, curr_mag ./ prev_mag, curr_energy ./ prev_energy, prev_x];
            featureV(end+1,:) = feature;
        end
    end
end
featureV(1,:) = [];
writematrix(featureV,'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\12-1-21_localization_dataset2.csv') 

%% split train/test set 12/2/21

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\';
takes = {'regular1', 'regular2', 'slow', 'stiffRknee', 'stiffLknee', 'insoleRweightL1', 'insoleRweightL2'};

first_stepV = zeros(1,16);
featureV = zeros(1,22);
for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie, intervention];
    load(string(filename))
    
    for i = 2:length(impacts)
        xcoord = coordinates(i,1)/1000;
        curr_mag = peak_mag(i,:);
        curr_energy = energy_squared(i,:);
        min_arr = min(arrival_idx(i,:));
        arrivals = arrival_idx(i,:) - min_arr;
        
        if walk_segments(i) ~= -1 % don't count prev step if it's a new walking segment
            prev_x = coordinates(i-1,1)/1000;
            prev_mag = peak_mag(i-1,:);
            prev_energy = energy_squared(i-1,:);
            
            % xy coord, arrival idx, peak mag, energy, xy prev coord
            feature = [xcoord, arrivals, curr_mag, curr_energy, curr_mag ./ prev_mag, curr_energy ./ prev_energy, prev_x];
            featureV(end+1,:) = feature;
        else
            first_stepV(end+1,:) = [xcoord, arrivals, curr_mag, curr_energy];
        end
    end
end
featureV(1,:) = [];
writematrix(featureV,'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\12-1-21_localization_dataset2.csv') 

%% for just the starts of segments 12/5/21

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\';
takes = {'regular1', 'regular2', 'slow', 'stiffRknee', 'stiffLknee', 'insoleRweightL1', 'insoleRweightL2'};

featureV = zeros(1,13);
for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie, intervention];
    load(string(filename))
    
    for i = 1:length(impacts)
        if walk_segments(i) == -1 % start of segment
            xcoord = coordinates(i,1)/1000;
            
            curr_mag = peak_mag(i,:);
            curr_energy = energy_squared(i,:);
            
            min_arr = min(arrival_idx(i,:));
            arrivals = arrival_idx(i,:) - min_arr;
            
            % xy coord, arrival idx, peak mag, energy, xy prev coord
            feature = [xcoord, arrivals, curr_mag, curr_energy];
            featureV(end+1,:) = feature;
        end
    end
end
featureV(1,:) = [];
writematrix(featureV,'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\12-5-21_localization_walkstarts.csv') 

%% taking predicted walk starts and creating second dataset 12/6/21

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\walkstarts_results.csv';
T = readtable(data_root_katie);
A = table2array(T);
A(1,:) = [];

idx = A(:,2);
realcoord = A(:,3);
predictcoord = A(:,4);
matrix = [idx,realcoord,predictcoord];
newmat = sortrows(matrix,1);

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\';
takes = {'regular1', 'regular2', 'slow', 'stiffRknee', 'stiffLknee', 'insoleRweightL1', 'insoleRweightL2'};

featureV = zeros(1,23);
start_count = 1;
% start = -1, train = 0, test = 1
for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie, intervention];
    load(string(filename))
    
    Nimpacts = length(impacts);
    test_i = round(Nimpacts*0.25); % splits train/test in 75/25%
    
    for i = 2:length(impacts)
        xcoord = coordinates(i,1)/1000;
        curr_mag = peak_mag(i,:);
        curr_energy = energy_squared(i,:);
        min_arr = min(arrival_idx(i,:));
        arrivals = arrival_idx(i,:) - min_arr;
        
        if walk_segments(i) ~= -1 % don't count prev step if it's a new walking segment
            if walk_segments(i-1) == -1
                prev_x = newmat(start_count,3);
                start_count = start_count + 1;
            else
                prev_x = coordinates(i-1,1)/1000;
            end
            prev_mag = peak_mag(i-1,:);
            prev_energy = energy_squared(i-1,:);
            
            % xy coord, arrival idx, peak mag, energy, xy prev coord
            if i > test_i*3 % && i <= test_i*3
                feature = [1, xcoord, arrivals, curr_mag, curr_energy, curr_mag ./ prev_mag, curr_energy ./ prev_energy, prev_x];
            else
                feature = [0, xcoord, arrivals, curr_mag, curr_energy, curr_mag ./ prev_mag, curr_energy ./ prev_energy, prev_x];
            end
            featureV(end+1,:) = feature;
        end
    end
end
featureV(1,:) = [];
writematrix(featureV,'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\12-6-21_trackinglocalization_4thquartile.csv') 

%% taking predicted locations of all points and making GRF dataset 12/7/21

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\walkstarts_results.csv';
T = readtable(data_root_katie);
walkstarts = table2array(T);
walkstarts(1,:) = [];
walkstarts(:,1) = [];

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\trackingloc_1stquart_results.csv';
T = readtable(data_root_katie);
firstquart = table2array(T);
firstquart(1,:) = [];
firstquart(:,1) = [];

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\trackingloc_2ndquart_results.csv';
T = readtable(data_root_katie);
secquart = table2array(T);
secquart(1,:) = [];
secquart(:,1) = [];

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\trackingloc_3rdquart_results.csv';
T = readtable(data_root_katie);
thirdquart = table2array(T);
thirdquart(1,:) = [];
thirdquart(:,1) = [];

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\trackingloc_4thquart_results.csv';
T = readtable(data_root_katie);
fourthquart = table2array(T);
fourthquart(1,:) = [];
fourthquart(:,1) = [];

walkdata = [firstquart;secquart;thirdquart;fourthquart]; % 1st col is index (python numbering!!), 2nd col is predicted xloc



% idx = A(:,2);
% realcoord = A(:,3);
% predictcoord = A(:,4);
% matrix = [idx,realcoord,predictcoord];
% newmat = sortrows(matrix,1);
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\';
takes = {'regular1', 'regular2', 'slow', 'stiffRknee', 'stiffLknee', 'insoleRweightL1', 'insoleRweightL2'};

featureV = zeros(1,16);
totalNcount = 1;
walkstartcount = 1;
for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie, intervention];
    load(string(filename))
    
    for i = 2:length(impacts)
        if walk_segments(i) == -1 % start of segment
            walkstart_idx = find(walkstarts(:,1)+1 == walkstartcount);
            xcoord = walkstarts(walkstart_idx,3);
            walkstartcount = walkstartcount + 1;
        else
            walkdata_idx = find(walkdata(:,1)+1 == totalNcount);
%             xcoord = walkdata(walkdata_idx,2); % finding predicted walk value
            xcoord = walkdata(i+1,2);
        end
        if ~isempty(xcoord)
            ycoord = coordinates(i,2);
            s1 = [-3.590,-3.343];
            s2 = [-3.580,2.61];
            s3 = [3.639,2.11];
            s4 = [3.650,-3.412];
            dist1 = sqrt( (xcoord-s1(1)).^2 + (ycoord-s1(2)).^2 );
            dist2 = sqrt( (xcoord-s2(1)).^2 + (ycoord-s2(2)).^2 );
            dist3 = sqrt( (xcoord-s3(1)).^2 + (ycoord-s3(2)).^2 );
            dist4 = sqrt( (xcoord-s4(1)).^2 + (ycoord-s4(2)).^2 );

            accX = acc_pks(i,1);
            accY = acc_pks(i,2);
            accZ = acc_pks(i,3);
            acc_sumsq = sqrt(accX^2 + accY^2 + accZ^2);

            % take #, acc, dist 1-4, pk mag 1-4, energy 1-4
            feature = [take,abs(accY),dist1,dist2,dist3,dist4,peak_mag(i,:),energy_squared(i,:),xcoord,coordinates(i,1)./1000]; %last two are just to check
            featureV(end+1,:) = feature;
        end
        totalNcount = totalNcount + 1;
    end
end
featureV(1,:) = [];
newV = featureV(:,1:14);
% doesn't work yet, the rmse was 2.8 instead of 0.6 (12/7/21)
writematrix(featureV,'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\12_7_21_grf_dataset_predictedloc.csv') 










