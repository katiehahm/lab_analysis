% step time variation
steptimesR = [];
timeR = [];
steptimesL = [];
timeL = [];

for i = 2:length(impacts)
    if impacts(i,4) == 0 & impacts(i-1,4) == 1
        steptimesR(end+1) = (impacts(i,1)-impacts(i-1,1))/296.3;
        timeR(end+1) = impacts(i,1)/296.3;
    elseif impacts(i,4) == 1 & impacts(i-1,4) == 0
        steptimesL(end+1) = (impacts(i,1)-impacts(i-1,1))/296.3;
        timeL(end+1) = impacts(i,1)/296.3;
    end
end

elim = find(steptimesR > 1);
steptimesR(elim) = [];
timeR(elim) = [];
elim = find(steptimesL > 1);
steptimesL(elim) = [];
timeL(elim) = [];

%%
figure;
plot(timeL,steptimesL,'b*')
hold on
plot(timeR,steptimesR,'r*')
legend('Left step time','Right step time')
xlabel('Time (s)')
ylabel('Step time (s)')
% title('Observed regular walking step time of each foot')
title('Observed braced walking step time of each foot')