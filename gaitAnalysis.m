function [A,B] = gaitAnalysis(onset_idx,Fs,impactN)
% 8/3/20
% detects the stride length of footsteps
stepA = [];
stepB = [];
prevT = 0;

for i = 1:impactN
    onsetavg = sum(onset_idx(i,:))/3;
    if i == 1
        prevT = onsetavg;
    elseif mod(i,2) == 0
        stepA(end+1) = onsetavg - prevT;
        prevT = onsetavg;
    else
        stepB(end+1) = onsetavg - prevT;
        prevT = onsetavg;
    end
end

if mod(impactN,2) == 1
    A = (sum(stepA)*2/(impactN-1))/Fs;
    B = (sum(stepB)*2/(impactN-1))/Fs;
else
    A = (sum(stepA)*2/impactN)/Fs;
    B = (sum(stepB)/(impactN/2 -1))/Fs;
end

figure;
hold on
plot(stepA./Fs, '*-')
plot(stepB./Fs, '*-')
xlabel('Step Number')
ylabel('Stride Time (s)')
legend('First leg', 'Second leg')
fprintf('First stride: %d ; Second stride: %d \n', A, B);

end

