%% Example - Mixture of von Mises Fisher distributions
clear;

% Seed the random number generator
rng(7);

%% Generate data
n1 = 100;
n2 = 200;
n3 = 50;
mu1 = [0 0 1]'; 
mu2 = [1 0 0]';
mu3 = [0 1 0]';
k1  = 15;
k2 = 20;
k3 = 37.5;

data = [ vmfrnd(k1, mu1, n1); vmfrnd(k2, mu2, n2); vmfrnd(k3, mu3, n3) ]; 

%% Fit data
% Run a mixture model with 1 class only
mm1 = snob(data, {'vmf',1:3}, 'fixedstructure', true, 'k', 1);

% Print out the details of the model
mm_Summary(mm1);

% Run a mixture model with 2 classes
mm2 = snob(data, {'vmf',1:3}, 'fixedstructure', true, 'k', 2);

% Print out the details of the model
mm_Summary(mm2);

% Run a mixture model with 3 classes
mm3 = snob(data, {'vmf',1:3}, 'fixedstructure', true, 'k', 3);

% Print out the details of the model
mm_Summary(mm3);

% Let snob determine the amount of classes automatically
% Start with 8 classes
mm4 = snob(data, {'vmf',1:3}, 'k', 8);
mm_Summary(mm4);

%% Plot the data on a unit sphere
model = mm4;
[x,y,z] = sphere(30);           % draw a sphere
mesh(x, y, z, 'FaceAlpha', 0.2, 'LineStyle', '--', 'EdgeColor', '#E5E5E5');
[~,cid] = max(model.r,[],2);      % determine which class each point belongs to

hold on;
nClasses = max(cid);
cmap = hsv(nClasses*2); 
cmap = cmap(randperm(nClasses*2), :);
cmap(1:3,:) = [1 0 0; 0 1 0; 0 0 1];    % first 3 colors are red/green/blue
% Note: selecting colours should really be done with something like:
% https://www.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap
for i = 1:nClasses
    
    % Plot all data points in class 'i'
    ix = (cid == i);
    plot3(data(ix,1), data(ix,2), data(ix,3), '.', 'MarkerSize', 6, ...
        'color', cmap(i,:), 'Marker', 'o', 'MarkerFaceColor', cmap(i,:));
    
    % Plot the mean vector of class 'i'
    mu = model.class{i}.model{1}.theta(2:end);
    line([0 mu(1)], [0 mu(2)], [0 mu(3)], 'LineWidth', 3, 'color', cmap(i,:));
end
