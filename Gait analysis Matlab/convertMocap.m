function [mocapT, mocapL, mocapR] = convertMocap(T, Mmocap)
% takes datatable from mocap data and converts it into arrays
% 8/30/21

A = table2array(T);
mocapT = A(:,2);
mocapL = A(:,Mmocap('Lx'):Mmocap('Lz'));
mocapR = A(:,Mmocap('Rx'):Mmocap('Rz'));

end

