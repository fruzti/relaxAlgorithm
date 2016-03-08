function [micFreqData, srcFreqData, theta, eta, ...
    micTimeData, srcTimeData] = genTstMicData(K, p, L, trueDOA, trueTOA,...
    P,fs,f0,beta)
% function [micFreqData, srcFreqData, theta, eta, ...
%     micTimeData, srcTimeData] = genTstMicData(K, p, L, trueDOA, trueTOA,...
%     P,fs,f0,beta)
% --------------------------------------------------------------------
% Generate a pseudorandom signal
% micFreqData : freq-domain signal at the microphones
% srcFreqData : freq-domain signal at the transmitter
% theta : DOAs
% eta : TOAs
% micTimeData : time-domain signal at the microphones
% srcTimeData : time-domain signal at the transmitter
% K :  number of microphones
% p : frequency-radius equivalence
% L : number of signals (optional)
% trueDOA : input DOAs (optional)
% trueTOA : input DOAs (optional)
% P : number of periods for the sinusoid (optional)
% fs : sampling frequency (optional)
% beta : attenuation (optional)

    if nargin < 3
        L = 1;
    end
    
    if nargin > 3
        flag = 1;
    end
    
    if nargin < 6
        P = 5;
        fs = 10e3;
    end
    
    if nargin < 8
        f0 = 5e3;
    end

    if nargin < 9
        beta = ones(L,1);
    end
    % sinosoid frequency
    w = 2*pi*f0;
    % number of samples
    samples = ceil(fs * P/f0);

    t = (0:samples)/fs;
    % original signal
%     srcTimeData = sin(w*t).';
    srcTimeData = randn(samples+1,1);
    N = length(srcTimeData); l = floor(N/2);
    srcFreqData = applyFFT(srcTimeData,N);

    % simulation parameters
    if ~flag
        theta = zeros(L,1); eta = zeros(L,1);
        theta(1) = 2*pi*rand;
        eta(1) = 2*pi*rand;
    else
        theta = trueDOA; eta = trueTOA;
    end
    micTimeData = beta(1)*genDelayData(srcTimeData,theta(1),eta(1),K,p);
    
    % microphone data
    for i = 2:L
        
        if ~flag
            % random theta and eta for each component
            theta(i) = 2*pi*rand;
            eta(i) = 2*pi*randn;
        end
            
        micTimeData = micTimeData + ...
                beta(i)*genDelayData(srcTimeData,theta(i),eta(i),K,p);
    end

    micFreqData = getFreqMicData(micTimeData, N, K);
    
end