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

% close all
% data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj';
% takes = {'normal1', 'slow', 'insole', 'weight', 'count', 'box', 'normal2'};
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\';
% takes = {'slow'};
takes = {'regular1', 'regular2', 'slow','stiffRknee', 'stiffLknee', 'insoleRweightL1', 'insoleRweightL2'};
% noise_thresh_arr = [0.002,0.00125,0.0012,0.0015];
noise_thresh_arr = [0.0002,0.00015,0.00015,0.00025];
noise_window = Fs*0.2; % take the 0.2s before the arrival of impact to get noise

% for person = 1:10
%     subj = int2str(person);
take = 7;
%     for take = 1:length(takes)
%         if person < 6 && take == 6
%             continue
%         elseif person == 2 && take == 2
%             continue
%         else
            intervention = char(takes(take));
%             filename = [data_root_katie, subj, '_', intervention];
            filename = [data_root_katie,intervention];
            load(string(filename))
            Fs = 12800;
            energy = zeros(length(arrival_idx),4);
            energy_envelope = zeros(length(arrival_idx),4);
            window_end_idx = zeros(length(arrival_idx),4);
            ups = zeros(length(arrival_idx),4);
            for j = 1:4
                for i = 1:length(arrival_idx)
%                     noise_thresh = noise_thresh_arr(j);
                    if i == length(arrival_idx)
                        window_end = length(filt_pcbD);
                    else
                        window_end = arrival_idx(i+1,j);
                    end
                    curr_i = arrival_idx(i,j);
                    window = filt_pcbD(curr_i:window_end,j);
                    noise_thresh = max(abs(filt_pcbD(curr_i - noise_window:curr_i-noise_window*0.1)))+0.0000005;
%                     noise_thresh = noise_thresh_arr(j);
                    if ~isempty(window)
%                         [up,lo] = envelope(abs(window),300,'rms');
                        [up,lo] = envelope(window,300,'peak');
                        indeces = find(up < noise_thresh);
                        while isempty(indeces) % keep raising threshold until find the signal end
                            noise_thresh = noise_thresh + 0.000001;
                            indeces = find(up < noise_thresh);
                        end
                        window = window(1:indeces(1));
                        window_end_idx(i,j) = indeces(1) + curr_i;
                        window_energy = sum(abs(window).^2);
                        energy_envelope(i,j) = sum(up(1:indeces(1)).^2);
                        energy(i,j) = window_energy;
                    end
                end
                
            end
            figure;
            for n = 1:4
            subplot(4,1,n)
            plot(filt_pcbD(:,n),'k')
            hold on
            plot(arrival_idx(:,n),0,'r.','MarkerSize',10)
            plot(window_end_idx(:,n),0,'c.','MarkerSize',10)
            plot(arrival_idx(:,n), energy(:,n),'b.','MarkerSize',10)
            plot(arrival_idx(:,n), energy_envelope(:,n),'m.','MarkerSize',10)
            end
%         end
        energy_squared = energy;
        save(filename,'energy_squared','energy_envelope','-append')
%     end
% end

% remember to delete any 0 elements from matrix!! 
% soemtiems arrival_idx(i+1) is less than arrival_idx(i) so left those as 0

%% making GRF feature vector 11/11/21

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj';
takes = {'normal1', 'slow', 'insole', 'weight', 'count', 'box', 'normal2'};
python_root = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Python Data\';

% subj #, fsr value, x coord, y coord, peak mag x4, energy x4
grf_features = zeros(1,14);

for person = 1:10
    subj = int2str(person);

    for take = 1:length(takes)
        if person < 6 && take == 6
            continue
        elseif person == 2 && take == 2
            continue
        else
            intervention = char(takes(take));
            filename = [data_root_katie, subj, '_', intervention, '_extract_straight_paths'];
            load(string(filename))
            Fs = 12800;
            
            for i = 1:length(impacts)
                if whichfoot(i) == 0
                    xcoord = coordinates(i,1);
                    ycoord = coordinates(i,2);
                    s1 = [-3.590,-3.343];
                    s2 = [-3.580,0.261];
                    s3 = [3.639,0.211];
                    s4 = [3.650,-3.412];
                    dist1 = sqrt( (xcoord-s1(1)).^2 + (ycoord-s1(2)).^2 );
                    dist2 = sqrt( (xcoord-s2(1)).^2 + (ycoord-s2(2)).^2 );
                    dist3 = sqrt( (xcoord-s3(1)).^2 + (ycoord-s3(2)).^2 );
                    dist4 = sqrt( (xcoord-s4(1)).^2 + (ycoord-s4(2)).^2 );
                    grf_features(end+1,:) = [person,impacts(i,3),dist1,dist2,dist3,dist4,peak_mag(i,1:4),envelope_energy(i,1:4)];
                end
            end
        end
    end
end

grf_features( ~any(grf_features,2), : ) = [];  % delete zero rows
writematrix(grf_features,'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Python Data\grf_features_envelope_energy_extract_straight_paths_left.csv') 


%% GRF feature vector changing coordinates to distance to sensor 11/15/21

filepath = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Python Data\grf_features_extract_straight_paths.csv';
T = readtable(filepath);

A = table2array(T);
xcoord = A(:,3); % in m
ycoord = A(:,4);

% change these values after confirming:
s1 = [-3.590,-3.343];
s2 = [-3.580,0.261];
s3 = [3.639,0.211];
s4 = [3.650,-3.412];
dist1 = sqrt( (xcoord-s1(1)).^2 + (ycoord-s1(2)).^2 );
dist2 = sqrt( (xcoord-s2(1)).^2 + (ycoord-s2(2)).^2 );
dist3 = sqrt( (xcoord-s3(1)).^2 + (ycoord-s3(2)).^2 );
dist4 = sqrt( (xcoord-s4(1)).^2 + (ycoord-s4(2)).^2 );

new_grf_features = [A(:,1:2),dist1,dist2,dist3,dist4,A(:,5:end)];
writematrix(new_grf_features,'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Python Data\distance_to_sensor_grf_features_extract_straight_paths.csv') 

%% asymmetry with GRF 11/16/21

clear all
close all
% data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj';
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\';
subj = '2'; % number of subject
% takes = {'normal1', 'insole', 'weight','count', 'normal2'};
takes = {'regular1', 'regular2', 'slow', 'stiffRknee', 'stiffLknee', 'insoleRweightL1', 'insoleRweightL2'};

all_params = zeros(length(takes),4);

accel_bar_data = zeros(length(takes),2);

for take = 1:length(takes)
    intervention = char(takes(take));
%     filename = [data_root_katie, subj, '_', intervention, '_extract_straight_paths'];
    filename = [data_root_katie, intervention];

    load(string(filename))
    Fs = 12800;
    
    left_mag = [];
    right_mag = [];
    
    left_acc = [];
    right_acc = [];
    for i = 1:length(whichfoot)
        if whichfoot(i) == 1 % right foot
            right_mag(end+1) = impacts(i,3);
            right_acc(end+1) = acc_pks(i,2);
        elseif whichfoot(i) == 0 % left foot
            left_mag(end+1) = impacts(i,3);
            left_acc(end+1) = acc_pks(i,2);
        end
    end
    
    all_params(take,:) = [mean(left_mag),mean(right_mag),std(left_mag),std(right_mag)];
    accel_bar_data(take,:) = [mean(left_acc),mean(right_acc)];

end
all_params
bar_data = zeros(length(takes)-1,2);
for i = 1:length(takes)-1
    normal_1 = all_params(1,1);
    normal_r = all_params(1,2);
    bar_data(i,:) = [normal_1 - all_params(i+1,1),normal_r - all_params(i+1,2)];
end
figure;
bar(accel_bar_data)
set(gca,'XTickLabel',takes(1:end));
% title(['Subj ', subj, ' differences in GRF from 1st regular walk to interventions'])
title('Average tibial acceleration of L and R legs')
legend('Left foot','Right foot')


%% envelope energy extraction 11/17/21

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
            filename = [data_root_katie, subj, '_', intervention,'_extract_straight_paths'];
            load(string(filename))
            Fs = 12800;
            envelope_energy = zeros(length(arrival_idx),4);
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
                        up = up(1:indeces(1));
                        window_energy = sum(abs(up));
                        envelope_energy(i,j) = window_energy; % energy might need to be squared sig but still ok?
                    end
                end
            end
            save(filename,'envelope_energy','-append')
        end
    end
end

%% making just normal insole GRF feature vector 11/17/21

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj';
takes = {'normal1', 'insole'};
python_root = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Python Data\';

% subj #, fsr value, x coord, y coord, peak mag x4, energy x4
grf_features = zeros(1,14);

subj = '1';

for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie, subj, '_', intervention, '_extract_straight_paths'];
    load(string(filename))
    Fs = 12800;

    for i = 1:length(impacts)
        if whichfoot(i) == 1
            xcoord = coordinates(i,1);
            ycoord = coordinates(i,2);
            s1 = [-3.590,-3.343];
            s2 = [-3.580,0.261];
            s3 = [3.639,0.211];
            s4 = [3.650,-3.412];
            dist1 = sqrt( (xcoord-s1(1)).^2 + (ycoord-s1(2)).^2 );
            dist2 = sqrt( (xcoord-s2(1)).^2 + (ycoord-s2(2)).^2 );
            dist3 = sqrt( (xcoord-s3(1)).^2 + (ycoord-s3(2)).^2 );
            dist4 = sqrt( (xcoord-s4(1)).^2 + (ycoord-s4(2)).^2 );
            grf_features(end+1,:) = [person,impacts(i,3),dist1,dist2,dist3,dist4,peak_mag(i,1:4),envelope_energy(i,1:4)];
        end
    end
end

grf_features( ~any(grf_features,2), : ) = [];  % delete zero rows
writematrix(grf_features,'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Python Data\grf_features_subj1_normal_insole_right.csv') 


%% GRF vector with distance, energy, and accel data 11/23/21

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\';
takes = {'regular1', 'regular2', 'slow', 'stiffRknee', 'stiffLknee', 'insoleRweightL1', 'insoleRweightL2'};

featureV = [0,0,0,0,0,0,0,0,0,0,0,0,0,0];
for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie, intervention];
    load(string(filename))
    
    for i = 1:length(impacts)
        xcoord = coordinates(i,1);
        ycoord = coordinates(i,2);
        s1 = [-3590,-3343];
        s2 = [-3580,261];
        s3 = [3639,211];
        s4 = [3650,-3412];
        dist1 = sqrt( (xcoord-s1(1)).^2 + (ycoord-s1(2)).^2 )/1000;
        dist2 = sqrt( (xcoord-s2(1)).^2 + (ycoord-s2(2)).^2 )/1000;
        dist3 = sqrt( (xcoord-s3(1)).^2 + (ycoord-s3(2)).^2 )/1000;
        dist4 = sqrt( (xcoord-s4(1)).^2 + (ycoord-s4(2)).^2 )/1000;
        
        accX = acc_pks(i,1);
        accY = acc_pks(i,2);
        accZ = acc_pks(i,3);
        acc_sumsq = sqrt(accX^2 + accY^2 + accZ^2);
        
        % take #, acc, dist 1-4, pk mag 1-4, energy 1-4
        feature = [take,abs(accY),dist1,dist2,dist3,dist4,peak_mag(i,:),energy_squared(i,:)];
        featureV(end+1,:) = feature;
    end
end
featureV(1,:) = [];
% writematrix(featureV,'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\grf_features_yAccel_energyenvelope.csv') 
    
%% generate train/test data to be 20% of each 11/24/21    
    
separations = [1,105,230,359,493,624,755,863];
trains = [];
tests = [];

for i = 2:length(separations)
    p = randperm(separations(i)-separations(i-1));
    p = p + separations(i-1)-1;
    twentyp = round( length(p)*0.8); % changed it to 10%
    trains = [trains,p(1:twentyp)];
    tests = [tests,p(twentyp+1:end)];
end

new_featureV = [0,0,0,0,0,0,0,0,0,0,0,0,0,0];
for i = 1:length(trains)
    new_featureV(end+1,:) = [0,featureV(trains(i),2:end)];
end
for i = 1:length(tests)
    new_featureV(end+1,:) = [featureV(tests(i),:)];
end

new_featureV(1,:) = [];
writematrix(new_featureV,'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\grf_features_yAccel_trainingset80p.csv') 
    
    
%% wavelet transform 11/28/21

close all
clear all
load('C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\insoleRweightL1');
Fs = 12800;

clip = filt_pcbD(arrival_idx(1)-10000:arrival_idx(10),1);
figure; plot(clip)

figure;
[wt,f,coi] = cwt(clip,Fs);
cwt(clip,Fs)

xrec = icwt(wt,f,[60 1000],'SignalMean',mean(clip));
figure; plot(xrec)