function [] = plotData(d,t,loc_names, title_str)
    len = length(loc_names);
    figure;
    hold on
    for i = 1:len
        subplot(len, 1, i)
        plot(t, d(:,i))
        title(append(title_str, loc_names(i)))
        ylabel('Volts')
        xlabel('Time')
    end
    
    figure;
    hold on
    for i = 1:len
        subplot(len, 1, i)
        plot(d(:,i))
        title(append(title_str, loc_names(i)))
        ylabel('Volts')
        xlabel('Index')
    end
end