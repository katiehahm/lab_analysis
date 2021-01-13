function [] = plotFreqSpectrum(d, Fs, loc_names)
    len = length(loc_names);
    title_str = 'frequency spectrum ';
    figure;
    hold on
    for i = 1:len
        subplot(len, 1, i)
        [X, F] = pwelch(d(:,i),[],[],[],Fs);
        plot(F, pow2db(X))
        title(append(title_str, loc_names(i)))
        ylabel('dB')
        xlabel('Frequency (samples/sec)')
    end
end