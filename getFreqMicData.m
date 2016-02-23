function micFreqData = getFreqMicData(micTimeData,N,K)
% function micFreqData = getFreqMicData(micTimeData,N,K)
% --------------------------------------------------------
% micFreqData : fft coefficients of micFreqData
% micTimeData : time-domain signal of microphones

    micFreqData = zeros(size(micTimeData));
    
    for k = 1:K
        micFreqData(:,k) = N*applyFFT(micTimeData(:,k),N); 
    end
    
end
