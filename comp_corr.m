function [c] = comp_corr(d1, seg1, d2, seg2)
% 6/22/20
% computes the cross correlation between two data vectors
% Input: d1/d2 = two vectors of data, seg1/seg2 = indeces to separate
% footfall segments of data

    n1 = length(seg1)/2; % number of segments1
    n2 = length(seg2)/2; % number of segments2
    combs = n1*n2; % total number of combinations
    c = zeros(1, combs);
    for i = 1:n1
        d1_min = seg1(i*2-1);
        d1_max = seg1(i*2);
        for j = 1:n2
            d2_min = seg2(j*2-1);
            d2_max = seg2(j*2);
            [corr, lags] = xcorr(d1(d1_min:d1_max), d2(d2_min:d2_max));
            figure;
            stem(lags, corr)
            c(((i-1)*n2)+j) = max(corr);
        end
    end

end

