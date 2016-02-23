function estBeta = nfUpdateBeta(micTimeData, srcTimeData,micPos,estX,...
    estY, N, K, L,fs)
% function estBeta = nfUpdateBeta(micTimeData, srcTimeData,micPos,estX,...
%     estY, N, K, L,fs)
% ------------------------------------------------------------------------
% estBeta : updated beta
% micTimeData : time-domain signal from receivers
% srcTimeData : time-domain signal from transmitter
% micPos : microphone positions
% estX : x-axis position of sources
% estY : y-axis position of sources
% N : number of samples
% K : number of microphones
% L : number of sources
% fs : sampling frequency

    S = zeros(N*K,L);
    
    for l = 1:L
        srcPos = [estX(l); estY(l)];
        tmp = nfGenDelayData(srcTimeData,micPos,srcPos,fs);
        S(:,l) = tmp(:);
    end

    estBeta = abs(S\micTimeData(:));

end