function [mocapT, mocapL, mocapR] = convertMocap(T, Mmocap)
% takes datatable from mocap data and converts it into arrays
% 8/30/21

A = table2array(T);
mocapT = A(:,2);
mocapLraw = A(:,Mmocap('Lx'):Mmocap('Lz'));
mocapRraw = A(:,Mmocap('Rx'):Mmocap('Rz'));

% added 11/4/21 to interpolate missing data
mocapL = fillmissing(mocapLraw,'linear');
mocapR = fillmissing(mocapRraw,'linear');

end

