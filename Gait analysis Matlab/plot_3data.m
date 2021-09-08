function [] = plot_3data(pcbD,pcbT,fsrD,fsrT,mocapR,mocapL,mocapT,Mfsr)
% plots pcb,fsr, and mocap data for visualization
% plots Rheel, Rtoe, Lheel, Ltoe, pcb(all), mocapR(y-axis), mocapL(y-axis)
% in a 7 x 1 subplot
% Mfsr is a dictionary mapping the fsr location to array index

figure;
hold on

subplot(7,1,1)
plot(fsrT,fsrD(:,Mfsr('Rheel')))
xlabel('Time (s)')
ylabel('Volts')
title('FSR Right heel')

subplot(7,1,2)
plot(fsrT,fsrD(:,Mfsr('Rtoe')))
xlabel('Time (s)')
ylabel('Volts')
title('FSR Right toe')

subplot(7,1,3)
plot(fsrT,fsrD(:,Mfsr('Lheel')))
xlabel('Time (s)')
ylabel('Volts')
title('FSR Left heel')

subplot(7,1,4)
plot(fsrT,fsrD(:,Mfsr('Ltoe')))
xlabel('Time (s)')
ylabel('Volts')
title('FSR Left toe')

subplot(7,1,5)
plot(pcbT,pcbD)
xlabel('Time (s)')
ylabel('Volts')
title('Accelerometer')
legend('S1','S2','S3','S4')

subplot(7,1,6)
plot(mocapT,mocapR(:,2))
xlabel('Time (s)')
ylabel('Height [mm]')
title('Mocap right foot')

subplot(7,1,7)
plot(mocapT,mocapL(:,2))
xlabel('Time (s)')
ylabel('Height [mm]')
title('Mocap left foot')

end

