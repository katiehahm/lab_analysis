% to extract energy from an intervention
% 11/29/21

% close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\';
takes = {'regular1', 'regular2', 'slow','stiffRknee', 'stiffLknee', 'insoleRweightL1', 'insoleRweightL2'};
take = 1;
intervention = char(takes(take));
filename = [data_root_katie,intervention];
load(string(filename))
Fs = 12800;
impactN = length(arrival_idx);

energy = zeros(impactN,4);
energy_idx = zeros(impactN,4); % idx where shaking ends

for j = 1:4
    for i = 1:impactN
        if i == length(arrival_idx)
            window_end = length(filt_pcbD);
        else
            window_end = arrival_idx(i+1,j);
        end
        curr_i = arrival_idx(i,j);
        window = filt_pcbD(curr_i:window_end,j);
        flip_window = flip(window);
        
        flipped_end_idx = aic_pick(window, 'to_peak');
        window_end_idx = length(window) - flipped_end_idx + 1;
        window = window(1:window_end_idx);
        if ~isempty(window)
            [up,lo] = envelope(abs(window),300,'peak'); % rms
            window_energy = sum(abs(window).^2);
            energy_idx(i,j) = curr_i + window_end_idx;
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
plot(energy_idx(:,n),0,'c.','MarkerSize',10)
plot(arrival_idx(:,n), energy(:,n),'b.','MarkerSize',10)
end

% energy_squared = energy;
% save(filename,'energy_squared','energy_envelope','-append')

%% find average accel per foot 11/29/21

clear all
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\';
takes = {'regular1', 'regular2', 'slow', 'stiffRknee', 'stiffLknee', 'insoleRweightL1', 'insoleRweightL2'};

all_params = zeros(length(takes),4);

accel_bar_data = zeros(length(takes),2);

take = 2;
intervention = char(takes(take));
filename = [data_root_katie, intervention];

load(string(filename))
Fs = 12800;

left_acc = [];
right_acc = [];

for i = 1:length(whichfoot)
    if whichfoot(i) == 1 % right foot
        right_acc(end+1) = acc_pks(i,2);
    elseif whichfoot(i) == 0 % left foot
        left_acc(end+1) = acc_pks(i,2);
    end
end
mean(right_acc)
mean(left_acc)