function [estBeta, S]= nfUpdateBeta2(micTimeData, srcTimeData,micPos,estX,...
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

%     lambda = 0.005;
    
%     cvx_begin quiet
%         variable estBeta(L)
%         minimize ( norm(S*estBeta - micTimeData(:),2) + lambda*norm(estBeta,1) )
%     cvx_end
%     estBeta = abs(S\micTimeData(:));
    estBeta = (S\micTimeData(:));

end