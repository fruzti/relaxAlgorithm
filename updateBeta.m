function estBeta = updateBeta(micTimeData, srcTimeData, estDOA,...
    estTOA, N, K, L, rho)
% estBeta : updated beta
% micTimeData : time-domain signal from receivers
% srcTimeData : time-domain signal from transmitter
% estDOA : estimated DOAs
% estTOA : estimated TOAs
% N : number of samples
% K : number of microphones
% L : number of sources
% rho : radius of UCA

    S = zeros(N*K,L);
    
    for l = 1:L
        tmp = genDelayData(srcTimeData,estDOA(l), estTOA(l), K, rho);
        S(:,l) = tmp(:);
    end

    estBeta = S\micTimeData(:);

end