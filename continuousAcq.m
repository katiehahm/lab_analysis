function continuousAcq(src, event)
%     [data, timestamps, ~] = read(src, src.ScansAvailableFcnCount, "OutputFormat", "Matrix");
    d = event.Data;
    t = event.TimeStamps;
    subplot(4,1,1)
    hold on
    plot(t, d(:,1), 'b')
    subplot(4,1,2)
    hold on
    plot(t, d(:,2), 'b')  
    subplot(4,1,3)
    hold on
    plot(t, d(:,3), 'b')
    subplot(4,1,4)
    hold on
    plot(t, d(:,4), 'b')
    
    changeint = floor(t/10) - floor( (t-0.5)/10 );
    if changeint > 0 & t > 10
        close all
        figure('Renderer', 'painters', 'Position', [900 10 1000 800])
        hold on
    end
% plot(event.TimeStamps, event.Data, 'b');
%     global datas = [ datas; event.Data];
%     global times = [ times, event.TimeStamps];

% 1/8/21 commented this out
    persistent tempData;
    global datas
    if(isempty(tempData))
        tempData=[];
    end
    tempData = [tempData;event.Data];
    datas = tempData;
    
    persistent tTemp;
    global times
    if(isempty(tTemp))
        tTemp = [];
    end
    tTemp=[tTemp;event.TimeStamps];
    times = tTemp;
    
%     data = event.Data;
%     time = event.TimeStamps;
%     save log.mat data time
%     fprintf(fid,'%f,%f\n',event.TimeStamps,event.Data) %     if ~ishandle(ButtonHandle)
%         disp('Data collection terminated');
%         src.stop()
%     end
end