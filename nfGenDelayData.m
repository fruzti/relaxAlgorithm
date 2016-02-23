function sDelayed = nfGenDelayData(srcTimeData,micPos, srcPos, fs)
% function sDelayed = genDelayData(srcTimeData,micPos, srcPos)
% ---------------------------------------------------
% srcTimeData : time-domain source signal
% micPos : microphone positions
% srcPos : position of the source
    
    % length of signal
    N = length(srcTimeData);
    % number of harmonics
    l = floor(N/2);
    % number of microphones
    K = size(micPos,2);
    % fft of original signal
    srcFreqData = applyFFT(srcTimeData,N);
    % microphone signals init 
    sDelayed = zeros(N,K);
    % filling array
    a = nfGen_a(micPos, srcPos, l, N, fs);
    for k = 1:K
        sDelayed(:,k) = applyIFFT(a(:,k).*srcFreqData,N);
    end
    
end