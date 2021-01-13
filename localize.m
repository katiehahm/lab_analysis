function [y,x] = localize(data)
% hilbert, AIC, findpeaks, triangulation.
% finds the [x, y] location of strikes
% 02/14/20

    % testing for git 1/13/21
    
    y = [];
    x = [];
    dLB = hilbert(data(:,LB));
    dLT = hilbert(data(:,LT));
    dRB = hilbert(data(:,RB));
    dRT = hilbert(data(:,RT));
    seconds = round(length(dLB)/1e5);
    pkLB = [];
    locLB = [];
    pkLT = [];
    locLT = [];
    pkRB = [];
    locRB = [];
    pkRT = [];
    locRT = [];
    
    % AIC on 1 sec windows
    for s = 0:(seconds-1)
        i_min = s*1e5+1;
        i_max = (s+1)*1e5;
        if s == seconds - 1
            i_max = length(dLB);
        else
            [pk,loc] = findpeaks(-1.*AIC_calc(dLB(i_min:i_max)));
            pkLB = pkLB + pk;
            locLB = locLB + loc;
            
            [pk,loc] = findpeaks(-1.*AIC_calc(dLT(i_min:i_max)));
            pkLT = pkLT + pk;
            locLT = locLT + loc;
            
            [pk,loc] = findpeaks(-1.*AIC_calc(dRB(i_min:i_max)));
            pkRB = pkRB + pk;
            locRB = locRB + loc;
            
            [pk,loc] = findpeaks(-1.*AIC_calc(dRT(i_min:i_max)));
            pkRT = pkRT + pk;
            locRT = locRT + loc;
        end
    end
    
    
end

