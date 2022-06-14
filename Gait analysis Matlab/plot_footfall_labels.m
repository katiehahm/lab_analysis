function [] = plot_footfall_labels(Y_predict,impacts,detected_starts,Fs_pcb)
% to be used in experiment4_processing3 for plotting estimated footfall
% labeling 6/10/22

figure;
plot(detected_starts(Y_predict == 1)./Fs_pcb,0.5,'ro')
hold on
plot(detected_starts(Y_predict == 2)./Fs_pcb,0.5,'rx')
real_o = find(floor(impacts(:,2)/10) == 1);
plot(impacts(real_o,1),0,'bo')
real_x = find(floor(impacts(:,2)/10) == 2);
plot(impacts(real_x,1),0,'bx')
ylim([-0.5,1])

end

