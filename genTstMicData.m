function [micFreqData, srcFreqData, theta, eta, ...
    micTimeData, srcTimeData] = genTstMicData(K, rho, L)
% Generate a mixture of 3 pure sines at 1kHz to test the algorithm 
% micTimeData : time-domain signal at the microphones
% srcTimeData : time-domain signal at the transmitter
% micFreqData : freq-domain signal at the microphones
% srcFreqData : freq-domain signal at the transmitter
% K :  number of microphones
% rho : radius of UCA
% L : number of signals

    if nargin < 3
        L = 1;
    end

    % sinosoid frequency
    f = 1e3;
    w = 2*pi*f;
    % number of samples
    P = 5;
    fs = 10e3;
    samples = ceil(fs * P/f);

    t = (0:samples)/fs;
    % original signal
    srcTimeData = sin(w*t).'; N = length(srcTimeData); l = floor(N/2);
    srcFreqData = applyFFT(srcTimeData,N);

    % simulation parameters
    theta = zeros(L,1); eta = zeros(L,1);
    theta(1) = 2*pi*rand;
    eta(1) = 2*pi*rand;
    micTimeData = genDelayData(srcTimeData,theta(1),eta(1),K,rho);
    
    % microphone data
    for i = 2:L
        % random theta and eta for each component
        theta(i) = 2*pi*rand;
        eta(i) = 2*pi*randn;
        micTimeData = micTimeData + ...
                genDelayData(srcTimeData,theta(i),eta(i),K,rho);
    end
    
    micFreqData = zeros(size(micTimeData));
    
    for k = 1:K
        micFreqData(:,k) = N*applyFFT(micTimeData(:,k),N); 
    end
    
end