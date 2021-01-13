function [c] = self_corr(d,segs)
% 6/22/20
% Input: d = data of one sensor, segs = indeces are start and end of each footfall
% Output: c = a vector containing the max cross-correlation between each
% footfall (ie. given 4 footfalls, c is length 6)

    n = length(segs)/2; % number of segments
    x = 1:n;
    M = combnk(x,2);
    combs = length(M(:,1)); % total number of conbinations
    c = zeros(1, combs);
    for i = 1:combs
       d0_min = segs(M(i,1)*2 -1);
       d0_max = segs(M(i,1)*2);
       d1_min = segs(M(i,2)*2-1);
       d1_max = segs(M(i,2)*2);
       [corr, lags] = xcorr(d(d0_min:d0_max), d(d1_min:d1_max));
       figure;
       stem(lags, corr)
       c(i) = max(corr);
    end

end

% initial manual testing
% load('data\barefoot_walking_06-19-2020_15-12.mat')
% barefoot_datas = datas;
% barefoot_times = times;
% load('data\sneakers_walking_06-19-2020_15-13.mat')
% sneakers_datas = datas;
% sneakers_times = times;
% 
% filt_barefoot_datas = lpf_data(barefoot_datas);
% filt_sneakers_datas = lpf_data(sneakers_datas);
% 
% b_A1 = filt_barefoot_datas(:,1);
% b_E1 = filt_barefoot_datas(:,2);
% b_A5 = filt_barefoot_datas(:,3);
% 
% s_A1 = filt_sneakers_datas(:,1);
% s_E1 = filt_sneakers_datas(:,2);
% s_A5 = filt_sneakers_datas(:,3);
% 
% ibA1_0 = 0;
% ibA1_1 = 0;
% 
% isA1_0 = 68520; %28830;
% isA1_1 = 72510; %33090;
% isE1_0 = 68660; %28830;
% isE1_1 = 71970; %33150;
% 
% % [c,lags] = xcorr(s_A1(isA1_0:isA1_1), s_E1(isE1_0:isE1_1));
% % [c,lags] = xcorr(s_A1(28830:33090), s_A1(68520:72510));
% [c,lags] = xcorr(s_A1(isA1_0:isA1_1), b_A1(43790:48640));
% stem(lags,c)