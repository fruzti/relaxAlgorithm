function [micTimeData, micFreqData, srcFreqData] = ...
    nfGenTstData(srcTimeData, micPos, srcPos, fs)
% function [micTimeData, micFreqData, srcFreqData] = ...
%     nfGenTstData(srcTimeData, micPos, srcPos, fs)
% --------------------------------------------------------
% micTimeData : time-domain microphones signal
% micFreqData : fft coefficients from microphones signal
% srcFreqData : fft coefficients from the source signal
% srcTimeData : time domain source signal
% micPos : position of the microphones
% srcPos : position of the source
% fs : sampling frequency

    K = size(micPos,1);
    L = size(srcPos,1);
    
    micTimeData = 0;
    
    for l = 1:L
        micTimeData = micTimeData + ...
            nfGenDelayData(srcTimeData, micPos(:,1:2)', srcPos(l,:)', fs);
    end
    
    N = length(srcTimeData);
    
    micFreqData = getFreqMicData(micTimeData, N, K);
    srcFreqData = applyFFT(srcTimeData, N);

end