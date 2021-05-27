function [mocapT, mocapL, mocapR] = convertMocap(T)
% takes datatable from mocap data and converts it into arrays
% arrays are clipped from after the first stomp
% 5/26/21

A = table2array(T);
mocapT = A(:,2);
Lloc = zeros(length(mocapT),3);
Rloc = zeros(length(mocapT),3);
for i = 1:length(mocapT)
    Lx = ( A(i,3)+ A(i,6) )/2;
    Ly = A(i,4);
    Lz = ( A(i,5)+ A(i,8) )/2;
    Rx = ( A(i,12)+ A(i,15) )/2;
    Ry = A(i,13);
    Rz = ( A(i,14) + A(i,17) )/2;
    Lloc(i,:) = [Lx,Ly,Lz];
    Rloc(i,:) = [Rx,Ry,Rz];
end

[~,i] = max(Rloc(:,2));
[~,mocapSlam] = min(Rloc(i:i+500,2)); % find index where foot fell and hit floor
mocapSlam = mocapSlam+i-1;
mocapSlamT = mocapT(mocapSlam);

mocapT = linspace(0,mocapT(end)-mocapT(mocapSlam),length(mocapT)-mocapSlam+1);
mocapL = Lloc(mocapSlam:end,:);
mocapR = Rloc(mocapSlam:end,:);

end

