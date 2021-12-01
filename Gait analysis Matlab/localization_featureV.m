% to make excel for localization 11/30/21

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\';
takes = {'regular1', 'regular2', 'slow', 'stiffRknee', 'stiffLknee', 'insoleRweightL1', 'insoleRweightL2'};

featureV = zeros(1,10);
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
            
            % xy coord, arrival idx, peak mag, energy, xy prev coord
            feature = [xcoord,curr_mag ./ prev_mag, curr_energy ./ prev_energy, prev_x];
            featureV(end+1,:) = feature;
        end
    end
end
featureV(1,:) = [];
writematrix(featureV,'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\12-1-21_localization_dataset.csv') 