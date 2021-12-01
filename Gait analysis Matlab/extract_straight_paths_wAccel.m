% 11/23/21
% relies on the user getting off the floor each time
% uses the amount of time passed between impacts to decide when the turns
% exist. Used for data collected 11/21/21

clear all
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\';
takes = {'regular1', 'regular2', 'stiffRknee', 'stiffLknee', 'insoleRweightL1', 'insoleRweightL2'};
% takes = {'slow'};
Fs = 12800;
turn_thresh = 1.3;

for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie,intervention];
    load(filename)

    % assume there is a long delay where the "off floor" peaks have been
    % filtered out
    % walk_segments: start of segment is -1, end of segment is 1, 0 in between
    Nimpacts = length(impacts(:,1));
    walk_segments = zeros(Nimpacts,1);
    walk_segments(1) = -1;
    for i = 2:Nimpacts-1
        time_btw = ( arrival_idx(i,1)-arrival_idx(i-1,1) )/Fs;
        if time_btw > turn_thresh % change in direction
            walk_segments(i) = -1;
            walk_segments(i-1) = 1;
        end
    end
    walk_segments(Nimpacts) = 1;
    
    % visualize
    figure;
    plot(filt_pcbD(:,1))
    hold on
    starts = find(walk_segments == -1);
    lasts = find(walk_segments == 1);
    plot(arrival_idx(starts,1),0,'mx','MarkerSize',15)
    plot(arrival_idx(lasts,1),0,'bx','MarkerSize',15)
    plot(arrival_idx(:,1),0,'r.','MarkerSize',10)
    
    save(filename,'walk_segments','-append')
end
