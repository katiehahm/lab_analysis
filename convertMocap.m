function [mocapT, mocapL, mocapR] = convertMocap(T, Mmocap)
% takes datatable from mocap data and converts it into arrays
% arrays are clipped from after the first stomp
% 5/26/21

A = table2array(T);
mocapT = A(:,2);
Lloc = zeros(length(mocapT),3);
Rloc = zeros(length(mocapT),3);
for i = 1:length(mocapT)
    Lxarr = [A(i,Mmocap('Lhx')), A(i,Mmocap('Lbx'))];
    Lzarr = [A(i,Mmocap('Lhz')), A(i,Mmocap('Lbz'))];
    Rxarr = [A(i,Mmocap('Rhx')), A(i,Mmocap('Rbx'))];
    Rzarr = [A(i,Mmocap('Rhz')), A(i,Mmocap('Rbz'))];
    
    Lx = mean(Lxarr, 'omitnan');
    if isnan(A(i,Mmocap('Lhy')))
        Ly = A(i,Mmocap('Lby'));
    else
        Ly = A(i,Mmocap('Lhy'));
    end
    Lz = mean(Lzarr, 'omitnan');
    Rx = mean(Rxarr, 'omitnan');
    if isnan(A(i,Mmocap('Rhy')))
        Ry = A(i,Mmocap('Rby'));
    else
        Ry = A(i,Mmocap('Rhy'));
    end
    Rz = mean(Rzarr, 'omitnan');
    Lloc(i,:) = [Lx,Ly,Lz];
    Rloc(i,:) = [Rx,Ry,Rz];
end

figure;
plot(mocapT,Lloc(:,1))
title('Left')
figure;
plot(mocapT,Rloc(:,1))
title('right')

[~,i] = max(Rloc(:,2));
[~,mocapSlam] = min(Rloc(i:i+500,2)); % find index where foot fell and hit floor
mocapSlam = mocapSlam+i-1;
mocapSlamT = mocapT(mocapSlam);

mocapT = linspace(0,mocapT(end)-mocapT(mocapSlam),length(mocapT)-mocapSlam+1);
mocapL = Lloc(mocapSlam:end,:);
mocapR = Rloc(mocapSlam:end,:);

end
