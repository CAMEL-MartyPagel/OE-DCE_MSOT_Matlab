function [Y f NFFT] = calc_fft(sig)

sigM = mean(sig);
if (sigM > 1)
    sig = sig - ones(size(sig,1),1)*sigM;
end;

Fs = 40000000;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = 2030;                     % Length of signal
t = (0:L-1)*T;                % Time vector

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(sig,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

