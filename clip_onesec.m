function [pcbD,pcbT,fsrT,fsrD,mocapT,mocapR,mocapL] = clip_onesec(pcbD,pcbT,fsrT,fsrD,mocapT,mocapR,mocapL)
% clips the first 1 sec of the final datasets to erase any effects of the
% big stomp

i = findTindex(1,pcbT);
pcbT = pcbT(i:end);
pcbD = pcbD(i:end,:);

i = findTindex(1,fsrT);
fsrT = fsrT(i:end);
fsrD = fsrD(:,i:end);
fsrD = fsrD.'; % transpose this because it's confusing

i = findTindex(1,mocapT);
mocapT = mocapT(i:end);
mocapR = mocapR(i:end,:);
mocapL = mocapL(i:end,:);

end

