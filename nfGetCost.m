function cost = nfGetCost(micTimeData, srcTimeData, micPos, estX,...
    estY, estBeta, L, fs)
% function cost = nfGetCost(micTimeData, srcTimeData, estX,...
%     estY, estBeta, L, K, fs)
% ------------------------------------------------------------------------
% cost : cost of error when micTimeData is approximate by L signals 
%       with parameters estDOA, estTOA and estBeta
% micTimeData : time-domain data received at the microphones
% srcTimeData : time-domain data from the transmitter
% micPos : microphones positions
% estX : estimated x-axis positions
% estY : estimated y-axis positions
% estBeta : estimated amplitude
% L : number of sources
% fs: sampling frequency
    
    yEst = 0;
    
    for l = 1 : L
        srcPos = [estX(l); estY(l)];
        yEst = yEst + estBeta(l) * nfGenDelayData(srcTimeData, micPos,...
            srcPos,fs);
    end
    
    cost = norm(micTimeData - yEst) ^ 2;


end