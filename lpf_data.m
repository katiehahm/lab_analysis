function [d] = lpf_data(datas)
%     beta = 3;
%     B = 50; % width of each band in Hz
%     L = 4; % length of each FIR filter
%     kw = kaiser(L, beta);
%     lpf = fir1(L-1, B/Fs, 'low', chebwin(35,30));

%     order = 34; % order of filter
%     cutfreq = 500*2/Fs; % desired cutoff freq * 2 / Fs
%     % Chebyshev window with 30 dB of ripple
%     lpf = fir1(order, cutfreq, 'low', chebwin(35,30));

    lpf = fir1(48, 0.3, 'low');
%     figure;
%     freqz(lpf,1)
    d = filter(lpf,1,datas);
end