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

%% Reduce Fs, Hilbert
reducedFs_data = downsample(datas,10); % keep every 10th datapoint
hilbert_data = imag(hilbert(reducedFs_data));

figure;
plot(hilbert_data(:,1));

noise = hilbert_data(65637:67172,1);
largeimpact = hilbert_data(67172:67737,1);
smallimpact = hilbert_data(73446:73999,1);

figure;
plot(noise)
title('Noise')
figure;
plot(largeimpact)
title('Large impact')
figure;
plot(smallimpact)
title('Small impact')

Fs = Fs/10;

y = noise;
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
title('Magnitude spectrum of noise in frequency');

y = largeimpact;
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
title('Magnitude spectrum of large impact in frequency');

y = smallimpact;
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
title('Magnitude spectrum of small impact in frequency');
% 
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

%% trying to cross correlate prev and after signals with each other 9/28
seg1 = hilbert_data(67001:67960,1);
seg2 = hilbert_data(67960:68687,1);
seg3 = hilbert_data(68687:69347,1);
seg4 = hilbert_data(69347:70117,1);

one2 = xcorr(seg1,seg2);
one3 = xcorr(seg1,seg3);
one4 = xcorr(seg1,seg4);

figure;
subplot(3,1,1)
plot(one2)
subplot(3,1,2)
plot(one3)
subplot(3,1,3)
plot(one4)

two1 = xcorr(seg2,seg1);
two3 = xcorr(seg2,seg3);
two4 = xcorr(seg2,seg4);

figure;
subplot(3,1,1)
plot(two1)
subplot(3,1,2)
plot(two3)
subplot(3,1,3)
plot(two4)

three1 = xcorr(seg3,seg1);
three2 = xcorr(seg3,seg2);
three4 = xcorr(seg3,seg4);
figure;
subplot(3,1,1)
plot(three1)
subplot(3,1,2)
plot(three2)
subplot(3,1,3)
plot(three4)

%% wavelet transform 9/28/21
cwt(hilbert_data(:,1),1)








