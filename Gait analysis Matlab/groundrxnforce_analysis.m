%% accumulating all forces 11/10/21
% just to visualize getting started

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj';
takes = {'normal1', 'slow', 'insole', 'weight', 'count', 'box', 'normal2'};
forces = zeros(1,5);
for j = 1:10
    subj = int2str(j);

    for take = 1:length(takes)
        if j < 6 && take == 6
            continue
        elseif j == 2 && take == 2
            continue
        else
            intervention = char(takes(take));
            filename = [data_root_katie, subj, '_', intervention];

            load(string(filename))
            Fs = 12800;

            % for extract_straight_paths:
        %     differences = [];
        %     for i = 1:(length(data)-1)
        %         if data(i+1,19) ~= -1 % edge doesn't say it's the start of seg
        %             curr = min(data(i+1,1:4));
        %             prev = min(data(i,1:4));
        %             differences(end+1) = (curr - prev)/Fs;
        %         end
        %     end

            % for original data
            for i = 1:length(impacts)
                forces(end+1,:) = [impacts(i,3),peak_mag(i,:)];
            end
        end
        
    end
end

% plotting heel impact forces in relation to pcb peak mag

figure;
plot(forces(:,1),max(forces(:,2:5),[],2),'.')
ylim([0 2.5])

figure;
plot(forces(:,1),'.')

%% energy extraction 11/10/21

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj';
takes = {'normal1', 'slow', 'insole', 'weight', 'count', 'box', 'normal2'};
noise_thresh_arr = [0.002,0.00125,0.0012,0.0015];


for person = 1:10
    subj = int2str(person);

    for take = 1:length(takes)
        if person < 6 && take == 6
            continue
        elseif person == 2 && take == 2
            continue
        else
            intervention = char(takes(take));
            filename = [data_root_katie, subj, '_', intervention];
            load(string(filename))
            Fs = 12800;
            energy = zeros(length(arrival_idx),4);
            for j = 1:4
                for i = 1:length(arrival_idx)
                    noise_thresh = noise_thresh_arr(j);
                    if i == length(arrival_idx)
                        window_end = length(filt_pcbD);
                    else
                        window_end = arrival_idx(i+1,j);
                    end
                    window = filt_pcbD(arrival_idx(i,j):window_end,j);
                    if ~isempty(window)
                        [up,lo] = envelope(abs(window),500,'rms');
                        indeces = find(up < noise_thresh);
                        while isempty(indeces) % keep raising threshold until find the signal end
                            noise_thresh = noise_thresh + 0.0005;
                            indeces = find(up < noise_thresh);
                        end
                        window = window(1:indeces(1));
                        window_energy = sum(abs(window));
                        energy(i,j) = window_energy; % energy might need to be squared sig but still ok?
                    end
                end
            end
        end
        save(filename,'energy','-append')
    end
end

% remember to delete any 0 elements from matrix!! 
% soemtiems arrival_idx(i+1) is less than arrival_idx(i) so left those as 0

%% making GRF feature vector 11/11/21

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj';
takes = {'normal1', 'slow', 'insole', 'weight', 'count', 'box', 'normal2'};
python_root = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Python Data\';

% subj #, fsr value, x coord, y coord, peak mag x4, energy x4
grf_features = zeros(1,12);

for person = 1:10
    subj = int2str(person);

    for take = 1:length(takes)
        if person < 6 && take == 6
            continue
        elseif person == 2 && take == 2
            continue
        else
            intervention = char(takes(take));
            filename = [data_root_katie, subj, '_', intervention];
            load(string(filename))
            Fs = 12800;
            
            for i = 1:length(impacts)
                grf_features(end+1,:) = [person,impacts(i,3),coordinates(i,1:2),peak_mag(i,1:4),energy(i,1:4)];
            end
        end
    end
end


grf_features( ~any(grf_features,2), : ) = [];  % delete zero rows
writematrix(grf_features,'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Python Data\grf_features.csv') 


%% GRF feature vector changing coordinates to distance to sensor 11/15/21

filepath = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Python Data\grf_features.csv';
T = readtable(filepath);

A = table2array(T);
xcoord = A(:,3); % in m
ycoord = A(:,4);

% change these values after confirming:
s1 = [-3590,-3343];
s2 = [-3580,261];
s3 = [3639,211];
s4 = [3650,-3412];
dist1 = sqrt( (xcoord-s1(1)).^2 + (ycoord-s1(2)).^2 );
dist2 = sqrt( (xcoord-s2(1)).^2 + (ycoord-s2(2)).^2 );
dist3 = sqrt( (xcoord-s3(1)).^2 + (ycoord-s3(2)).^2 );
dist4 = sqrt( (xcoord-s4(1)).^2 + (ycoord-s4(2)).^2 );

new_grf_features = [A(:,1:2),dist1,dist2,dist3,dist4,A(:,5:end)];
writematrix(new_grf_features,'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Python Data\distance_to_sensor_grf_features.csv') 




