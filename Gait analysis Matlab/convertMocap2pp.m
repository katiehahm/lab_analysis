function [mocapT, mocapL1, mocapR1, mocapL2, mocapR2] = convertMocap2pp(T, Mmocap)
% takes datatable from mocap data and converts it into arrays with two
% people
% 2/8/22

A = table2array(T);
mocapT = A(:,2);
mocapLraw1 = A(:,Mmocap('Lx1'):Mmocap('Lz1'));
mocapRraw1 = A(:,Mmocap('Rx1'):Mmocap('Rz1'));
mocapLraw2 = A(:,Mmocap('Lx2'):Mmocap('Lz2'));
mocapRraw2 = A(:,Mmocap('Rx2'):Mmocap('Rz2'));

% added 11/4/21 to interpolate missing data
mocapL1 = fillmissing(mocapLraw1,'linear');
mocapR1 = fillmissing(mocapRraw1,'linear');
mocapL2 = fillmissing(mocapLraw2,'linear');
mocapR2 = fillmissing(mocapRraw2,'linear');

end

