function [xest,B,MSE] = wienerFilt(x,y,N)
%
% Wiener filter based on Wiener-Hopf equations
%   This function takes as inputs a noisy signal, x, and a reference signal, y,
%   in order to compute a N-order linear filter that provides an estimate of y
%   from x
%  
% INPUTS
% x = noisy signal
% y = reference signalsinit
% N = filter order
%
% OUTPUTS
% xest = estimated signal
% b = Wiener filter coefficents
% MSE = mean squared error
%
% M. Buzzoni
% May 2019
% -------------
% Rev. Feb. 2020: the function can be performed by using column or row
% vectors as inputs
X = 1/N .* fft(x(1:N));
Y = 1/N .* fft(y(1:N));
X = X(:);
Y = Y(:);
Rxx = N .* real(ifft(X .* conj(X))); % Autocorrelation function
Rxy = N .* real(ifft(X .* conj(Y))); % Crosscorrelation function
Rxx = toeplitz(Rxx);
Rxy = Rxy';
B = Rxy / Rxx; B = B(:); % Wiener-Hopf eq. B = inv(Rxx) Rxy
xest = fftfilt(B,x);
xest = xest(N+1:end); % cut first N samples due to distorsion during filtering operation
MSE = mean(y(N+1:end) - xest) .^2; % mean squared error
