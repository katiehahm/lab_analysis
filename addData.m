function [data] = addData(prevData, diff, mag_diff)
% 7/28/20
% forms a data bank for machine learning
d = [diff,mag_diff];
data = [prevData;d];

end

