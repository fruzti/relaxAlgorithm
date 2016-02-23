function estBeta = updateBeta(micTimeData, srcTimeData, estDOA,...
    estTOA, N, K, L, p)
% function estBeta = updateBeta(micTimeData, srcTimeData, estDOA,...
%     estTOA, N, K, L, p)
% ------------------------------------------------------------------------
% estBeta : updated beta
% micTimeData : time-domain signal from receivers
% srcTimeData : time-domain signal from transmitter
% estDOA : estimated DOAs
% estTOA : estimated TOAs
% N : number of samples
% K : number of microphones
% L : number of sources
% p : frequency-radius equivalence

    S = zeros(N*K,L);
    
    for l = 1:L
        tmp = genDelayData(srcTimeData,estDOA(l), estTOA(l), K, p);
        S(:,l) = tmp(:);
    end

    estBeta = abs(S\micTimeData(:));

end