function [clipped_time, clipped_data] = clip_to_fsr(time, data)
% uses triggers in fsr to clip both fsr and pcb datasets
% 5/26/21

invtrig = -1 .* data(:,5);
edges = [0,0,0,0];
for n = 1:4
    [~,i] = max(invtrig);
    edges(n) = i;
    invtrig(i - 10000:i+10000) = -5;
end
edges = sort(edges);
startI = edges(1);
endI = edges(3);

clipped_data = data(startI:endI,1:4);
clipped_time = linspace(0, time(endI)-time(startI), endI-startI+1);

end