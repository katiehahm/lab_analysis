load('data/data_10-02-2020_19-03.mat')
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

% [peaks1, peaks2, peaks3, idx1, idx2, idx3] = TDOA4new(clean_data) %try creating TDOA4
% d = clean_data;
% sn = length(d(1,:)); % number of sensors
% figure;
%     for i = 1:sn
%         hold on
%         subplot(sn,1,i)
%         title(append('Filtered data at ', loc_names(1)))
%         xlabel('Impact number')
%         ylabel('Volts (V)')
%         hold on
%         plot(d(:,i))
% %         plot(onset_idx(i,:),0,'bx')
% %         plot(idx(i,:),peaks(i,:),'ro')
%         if i ==1
%             hold on
%             plot(idx1,peaks1,'ro')
%         end
%         if i ==2
%             hold on
%             plot(idx2,peaks2,'ro')
%         end
%         if i ==3
%             hold on
%             plot(idx3,peaks3,'ro')
%         end
%     end


%     for i = 1:sn
%         hold on
%         subplot(sn,1,i)
%         title(append('Filtered data at ', loc_names(1)))
%         xlabel('Impact number')
%         ylabel('Volts (V)')
%         hold on
%         plot(d(i,:))
% %         plot(onset_idx(i,:),0,'bx')
% %         plot(idx(i,:),peaks(i,:),'ro')
%         if i ==1
%             hold on
%             plot(idx1,peaks1,'ro')
%         end
%         if i ==2
%             hold on
%             plot(idx2,peaks2,'ro')
%         end
%         if i ==3
%             hold on
%             plot(idx3,peaks3,'ro')
%         end
%         
%     end

impactN = 16; % variable input value

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
[cadencetotal,avgstridelengths,stddevstridelength] = cadence(peak_idx,impactN,Fs)


avgcadence = [mean([1.5715 1.6656 1.4359]) mean([1.5956 1.5950 1.5989]) mean([1.8839 1.8680 1.8762]) mean([1.7490 1.7202 1.7332]) mean([1.5211 1.5224 1.5229]) mean([1.4414 1.4378 1.4320]) mean([1.3811 1.3946 1.0988]) mean([1.3514 1.3486 1.0414])]
avgstridelength = [0.6940 0.6711 0.6397 0.6920 0.7117 0.7539 0.8327 0.8879]
limpingc = avgcadence(6:8)
limpings=avgstridelength(6:8)
carolync = avgcadence(1:2)
carolyns= avgstridelength(1:2)
figure
scatter(avgcadence, avgstridelength)
hold on
scatter(limpingc,limpings,'r')
scatter(carolync,carolyns,'g')
title('Stride Length vs. Cadence')
xlabel('Average Cadence')
ylabel('Average Stride Length')
legend('Katie walking','Katie limping', 'Carolyn walking')









