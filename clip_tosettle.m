function [pcbD,pcbT,fsrT,fsrD,mocapT,mocapR,mocapL] = clip_tosettle(timeseg,pcbD,pcbT,fsrT,fsrD,mocapT,mocapR,mocapL)
% clips the first 1 sec of the final datasets to erase any effects of the
% big stomp

i = findTindex(timeseg,pcbT);
pcbT = pcbT(i:end);
pcbD = pcbD(i:end,:);

i = findTindex(timeseg,fsrT);
fsrT = fsrT(i:end);
fsrD = fsrD(:,i:end);
fsrD = fsrD.'; % transpose this because it's confusing

i = findTindex(timeseg,mocapT);
mocapT = mocapT(i:end);
mocapR = mocapR(i:end,:);
mocapL = mocapL(i:end,:);

end

