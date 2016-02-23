function sDelayed = genDelayData(s,theta,eta,K,p)
% function sDelayed = genDelayData(s,theta,eta,K,p)
% ---------------------------------------------------
% sDelay : delayed time signal
% theta : direction of arrival
% eta : angle of range
% K : number of microphones
% p : frequency-radius equivalence
    
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
        dk = gen_d(theta,l,K,k,p);
        sDelayed(:,k) = applyIFFT(a.*dk.*f,N);
    end
end