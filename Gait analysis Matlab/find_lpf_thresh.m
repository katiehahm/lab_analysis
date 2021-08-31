data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\';
datestr = '06-10-2021_14-17-29';
load([data_root_katie, 'ImmersionFloorSensors\data_', datestr])
plot(times, datas(:,1))
plot(datas(:,1))
data = datas(664101:end,:);
time = times(664101:end);
close all
plot(time, data(:,1))
plot(data(:,1))
footstep1 = data(6657:13476,1);
noise1 = data(13476:16293);
Fs = 12800;

figure;
plot(footstep1)
title('Footstep time signal')

figure;
plot(noise1)
title('Noise time signal')

y = footstep1;
N = length(y);              % Length of vector y, number of samples
Y = fft(y,N);               % Fourier transform of y
F = ((0:1/N:1-1/N)*Fs);     % Frequency vector
w = 2*pi*F;                 % Angular frequency vector
magnitudeY = abs(Y);        % Magnitude of the FFT
phaseY = unwrap(angle(Y));  % Phase of the FFT
figure;
% subplot(2,1,1);
plot(F, magnitudeY);
footstep1_mag = magnitudeY;
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Frequency, Hz');
ylabel('Magnitude, dB');
title('Magnitude spectrum of footstep in frequency');
% subplot(2,1,2);
% plot(w, magnitudeY);
% grid on;
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
% xlabel('Angular Frequency, rad/sample');
% ylabel('Magnitude, dB');
% title('Magnitude spectrum of sound wave in angular frequency');
% figure;
% subplot(2,1,1);
% plot (F, phaseY);
% grid on;
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
% xlabel('Frequency, Hz');
% ylabel('Phase angle, radian');
% title('Phase spectrum of sound wave in frequency');
% subplot(2,1,2);
% plot (w, phaseY);
% grid on;
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
% xlabel('Angular Frequency, rad/sample');
% ylabel('Phase angle, radian');
% title('Phase spectrum of sound wave in angular frequency');

y = noise1;
N = length(y);              % Length of vector y, number of samples
Y = fft(y,N);               % Fourier transform of y
F = ((0:1/N:1-1/N)*Fs);     % Frequency vector
w = 2*pi*F;                 % Angular frequency vector
magnitudeY = abs(Y);        % Magnitude of the FFT
phaseY = unwrap(angle(Y));  % Phase of the FFT
figure;
% subplot(2,1,1);
plot(F, magnitudeY);
noise1mag = magnitudeY;
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Frequency, Hz');
ylabel('Magnitude, dB');
title('Magnitude spectrum of noise in frequency');
% subplot(2,1,2);
% plot(w, magnitudeY);
% grid on;
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
% xlabel('Angular Frequency, rad/sample');
% ylabel('Magnitude, dB');
% title('Magnitude spectrum of sound wave in angular frequency');
% figure;
% subplot(2,1,1);
% plot (F, phaseY);
% grid on;
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
% xlabel('Frequency, Hz');
% ylabel('Phase angle, radian');
% title('Phase spectrum of sound wave in frequency');
% subplot(2,1,2);
% plot (w, phaseY);
% grid on;
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
% xlabel('Angular Frequency, rad/sample');
% ylabel('Phase angle, radian');
% title('Phase spectrum of sound wave in angular frequency');

footstep2 = data(513653:516581,1);
noise2 = data(516581:520714,1);
figure; plot(footstep2)
title('Footstep 2')
figure; plot(noise2)
title('Noise 2')

y = footstep2;
N = length(y);              % Length of vector y, number of samples
Y = fft(y,N);               % Fourier transform of y
F = ((0:1/N:1-1/N)*Fs);     % Frequency vector
w = 2*pi*F;                 % Angular frequency vector
magnitudeY = abs(Y);        % Magnitude of the FFT
phaseY = unwrap(angle(Y));  % Phase of the FFT
figure;
footstep2mag = magnitudeY;
plot(F, magnitudeY);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Frequency, Hz');
ylabel('Magnitude, dB');
title('Magnitude spectrum of footstep2 in frequency');

y = noise2;
N = length(y);              % Length of vector y, number of samples
Y = fft(y,N);               % Fourier transform of y
F = ((0:1/N:1-1/N)*Fs);     % Frequency vector
w = 2*pi*F;                 % Angular frequency vector
magnitudeY = abs(Y);        % Magnitude of the FFT
phaseY = unwrap(angle(Y));  % Phase of the FFT
figure;
plot(F, magnitudeY);
noise2mag = magnitudeY;
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Frequency, Hz');
ylabel('Magnitude, dB');
title('Magnitude spectrum of noise2 in frequency');