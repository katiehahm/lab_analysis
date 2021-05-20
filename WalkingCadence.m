function [cadence, stan, vel] = WalkingCadence(datas)
    time_between = [];
    stride = 70;
    Fs = 12800;
    for n = 1:length(datas) 
        data = datas{n};
        [peaks1, peaks2, peaks3, idx1, idx2, idx3] = TDOA4new(data);
        num_steps = length(peaks1);
        for i = 1:num_steps-1
            val = (idx1(i+1) - idx1(i))/Fs;
            time_between = [time_between, val]; 
        end 
    end
    div = 1./time_between;
    mult = 60*div;
    cadence = mean(mult);
    stan = std(mult);
    vel = cadence/60*stride*1e-2;
end