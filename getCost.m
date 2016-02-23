function cost = getCost(micTimeData, srcTimeData, estDOA,...
    estTOA, estBeta, L, K, p)
% function cost = getCost(micTimeData, srcTimeData, estDOA,...
%     estTOA, estBeta, L, K, p)
% ------------------------------------------------------------------------
% cost : cost of error when micTimeData is approximate by L signals 
%       with parameters estDOA, estTOA and estBeta
% micTimeData : time-domain data received at the microphones
% srcTimeData : time-domain data from the transmitter
% estDOA : estimated DOAs
% estTOA : estimated TOAs
% estBeta : estimated amplitude
% L : number of sources
% K : number of microphones
% p : frequency-radius equivalence
    
    yEst = 0;
    
    for l = 1 : L
        yEst = yEst + estBeta(l) * genDelayData(srcTimeData,estDOA(l),...
            estTOA(l),K, p);
    end
    
    cost = norm(micTimeData - yEst) ^ 2;


end