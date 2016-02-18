function y = applyFFT(x,N)
% y : centered fft of x
% x : time-domain signal
% N : length of signal

    if nargin < 2
        N = length(x);
    end
    
    y = 1/N * fftshift(fft(x));

end