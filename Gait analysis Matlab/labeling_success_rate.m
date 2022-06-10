function [] = labeling_success_rate(impacts, detected_starts, Y_predict, Fs_pcb, Fs_fsr, rate)
% 6/8/22
% used to calculate success rate of labeling from detected impacts
% used in decisionTree_footfall_classifier.m

detected_times = detected_starts./Fs_pcb;
success_count = 0;
for i = 1:length(impacts)
    real_time = impacts(i,1)/Fs_fsr;
    near_idx = find(abs(detected_times - real_time) < rate);
    real_label = impacts(i,4);
    if ~isempty(near_idx)
        if real_label == 1 | real_label == 2
            for j = 1:length(near_idx)
                if Y_predict(near_idx(j)) == 1
                    success_count = success_count + 1;
                end
            end
        else
            for j = 1:length(near_idx)
                if Y_predict(near_idx(j)) == 2
                    success_count = success_count + 1;
                end
            end
        end
    end
end
success_rate = success_count / length(impacts)
