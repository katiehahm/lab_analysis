% 9/9/21, 11/16/21
% uses change in direction in x-dir to determine walking episode separation
% plots the original points and overlays with x with the deleted points
clear all
close all
% data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj';
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\';
% subj = '10'; % number of subject
takes = {'regular1', 'regular2', 'slow', 'stiffRknee', 'stiffLknee', 'insoleRweightL1', 'insoleRweightL2'};
% takes = {'normal1', 'slow', 'insole', 'weight', 'count', 'box','normal2'};
deleted_pts = cell(length(takes),1);

for take = 1:length(takes)
    intervention = char(takes(take));
%     filename = [data_root_katie, subj, '_', intervention];
    filename = [data_root_katie,intervention];
    load(filename)

    % 10/20/21 
    % assumes that subj always walks in the long direction, so 
    % any change in x direction is a turn
    curr_x = sign( coordinates(2,1)-coordinates(1,1) );
    for i = 2:(length(coordinates)-1)
        new_x = sign( coordinates(i+1,1)-coordinates(i,1) );
        if curr_x ~= new_x % change in direction
            curr_x = new_x;
            deleted_pts{take}(end+1) = i;
        end
    end
    % compare these two figures with the two at the end of this section
    % the turning points should've been eliminated
    figure; 
    plot(deleted_pts{take},coordinates(deleted_pts{take},1),'bx')
    hold on
    plot(coordinates(:,1),'r.')
end
%% based on above plots, add extra points
% 11/16/21

deleted_pts{1} = sort([deleted_pts{1}, 114]);
deleted_pts{2} = sort([deleted_pts{2} ]);
deleted_pts{3} = sort([deleted_pts{3} , 62,63,64,65,101,126,116,117,115,139,159,158,160,161]);
deleted_pts{4} = sort([deleted_pts{4} ,1,148]);
deleted_pts{5} = sort([deleted_pts{5}, 26,28,138,153]);
deleted_pts{6} = sort([deleted_pts{6}, 25,52,79,66,38,94,93]);
deleted_pts{7} = sort([deleted_pts{7}, 147]);

% omit_pts{1} = [];
% omit_pts{2} = [];
% omit_pts{3} = [];
% omit_pts{4} = [];
% omit_pts{5} = [];
% % omit_pts{6} = [];


for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie, subj, '_', intervention];
    load(filename)
    
    figure; 
    plot(deleted_pts{take},coordinates(deleted_pts{take},1),'bx')
    hold on
    plot(coordinates(:,1),'r.')
    title(intervention)
end

%% then save the new data
% 11/16/21

for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie, subj, '_', intervention];
    load(filename)
    
    delete_idx = deleted_pts{take};
    
    
    % start of segment is -1, end of segment is 1, 0 in between
    walk_episodes = zeros(length(arrival_idx),1);
    if delete_idx(1) ~= 1
        walk_episodes(1) = -1;
        walk_episodes(delete_idx(1)-1) = 1;
    end
    for i = 2:length(delete_idx)
        % assume if more than one is deleted, they are always consecutive
        if delete_idx(i)-delete_idx(i-1) > 1 % not consecutive
            walk_episodes(delete_idx(i)-1) = 1;
            walk_episodes(delete_idx(i-1)+1) = -1;
        end
    end
    if delete_idx(end) < length(walk_episodes)
        walk_episodes(end) = 1;
        walk_episodes(delete_idx(end)+1) = -1;
    end
    arrival_idx(delete_idx,:) = [];
    coordinates(delete_idx,:) = [];
    energy(delete_idx,:) = [];
    impacts(delete_idx,:) = [];
    peak_idx(delete_idx,:) = [];
    peak_mag(delete_idx,:) = [];
    whichfoot(delete_idx,:) = [];
    walk_episodes(delete_idx,:) = [];
    
    figure;
    plot(coordinates(:,1),'k.')
    title(intervention)
    hold on
    plot(find(walk_episodes==1),coordinates(find(walk_episodes==1),1),'bx')
    plot(find(walk_episodes==-1),coordinates(find(walk_episodes==-1),1),'rx')
    
    filename = [filename, '_extract_straight_paths'];
    save(filename,'arrival_idx','coordinates','energy','impacts','peak_idx','peak_mag','whichfoot','walk_episodes')
end

%% forgot to add some variables to extract_straight_paths 11/17/21

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj';
takes = {'normal1', 'slow', 'insole', 'weight', 'count', 'box', 'normal2'};

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
            filename2 = [data_root_katie, subj, '_', intervention,'_extract_straight_paths'];
            load(string(filename))
            
            save(filename2,'filt_pcbD','fsrData','fsrTime','mocapL','mocapR','mocapT','pcbTime','-append')
        end
    end
end