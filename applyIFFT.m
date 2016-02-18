function x = applyIFFT(y,N)
% x : time-domain signal
% y : centered fft of x
% N : length of signal

    x = ifft(ifftshift(y));
    
    if nargin < 2
        N = length(x);
    end
    
    x = N * x;
end