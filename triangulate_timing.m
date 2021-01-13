function [mag_diff, diff] = triangulate_timing(idxs,vals, loc_names,Fs,plot)
% Performs triangulation on the given index and corresponding values
% Returns mag_diff: a matrix of the differences in magnitude
% Returns diff: a matrix of the differences in timing

if nargin < 5
    plot = true;
end

npeaks = length(idxs);
s = length(idxs(1,:));
diff = zeros(npeaks, s);
mag_diff = zeros(npeaks, s);

for i = 1:npeaks
    min_idx = idxs(i,1);
    min_mag = vals(i,1);
    for j = 2:s
        min_idx = min([min_idx, idxs(i,j)]);
        min_mag = min([min_mag, vals(i,j)]);
    end
    for j = 1:s
        diff(i,j) = idxs(i,j) - min_idx;
        mag_diff(i,j) = vals(i,j) - min_mag;
    end
end

% convert number of samples to time
diff = diff.*1000./Fs;

if (plot)
    figure;
    bar(diff)
    legend(loc_names)
    title('Triangulation: time delays of each sensor from peaks')
    xlabel('Impact #')
    ylabel('time(ms)')

    figure;
    bar(mag_diff)
    legend(loc_names)
    title('Triangulation: magnitude differences of each sensor at peak')
    xlabel('Impact #')
    ylabel('Volts')
end

end

