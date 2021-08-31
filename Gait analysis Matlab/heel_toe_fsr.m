function [RheelT, RtoeT, LheelT, LtoeT] = heel_toe_fsr(impactsT, fsrT, fsrD, Mfsr)
    % impacts are time stamps of the peaks of both L and R heels
    % ordered by increasing time
    
    heeltotoe = 0.1; % ########## change #######################
    
    N = length(impactsT);
    RheelT = zeros(1,N);
    RtoeT = zeros(1,N);
    LheelT = zeros(1,N);
    LtoeT = zeros(1,N);
    
    for i = 1:length(impactsT)
        % heel strike time is when going back in time heel mag is zero again
        
        % toe off time is when after time heeltotoe, it's 0 again
    end
end

