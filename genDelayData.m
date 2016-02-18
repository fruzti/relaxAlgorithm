function sDelayed = genDelayData(s,theta,eta,K,rho)
% sDelay : delayed time signal
% theta : direction of arrival
% eta : angle of range
% K : number of microphones
% rho : radius of array
    
    % length of signal
    N = length(s);
    % number of harmonics
    l = floor(N/2);
    
    % fft of original signal
    a = applyFFT(s,N);
    % range delay
    f = gen_f(eta,l);
    % microphone signals init 
    sDelayed = zeros(N,K);
    % filling array
    for k = 1:K
        dk = gen_d(theta,l,rho,K,k);
        sDelayed(:,k) = applyIFFT(a.*dk.*f,N);
    end
end