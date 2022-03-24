% used to detect footfalls
% 3/24/22

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Praneeth experiment 3_11_22\ProcessedData\both_towards1';
load(string(data_root_katie))


segments = [0, 18, 38, 59, 79, 101, 122, 143, 164, 185, 206, 225];
t = linspace(0,162199/Fs_pcb,162199);
figure;
impacttimes = impacts(impactidx,1)/Fs_fsr;
valididx = find(f < 600 & f > 100);
plot(impacttimes-pcbTime(113830),150,'rx','MarkerSize',8)
hold on
for i = 1:6
    pcbclip = pcbData(113830:276028,i);
    [wt,f] = cwt(pcbclip, Fs_pcb);
    hold on
    contour(t,f(valididx),abs(wt(valididx,:)))
%     set(gca, 'YScale', 'log')
end

