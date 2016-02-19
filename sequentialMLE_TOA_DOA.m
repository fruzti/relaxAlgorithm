function [estDOA, estTOA, estBeta] = sequentialMLE_TOA_DOA(micTimeData,...
    srcTimeData, srcFreqData, K, rho,L,N)
% estDOA : estimated DOAs
% estTOA : estimated TOAs
% estBeta : estimated Betas
% micTimeData : time-domain signal of the receivers
% srcTimeData : time-domain signal from the transmitter
% srcFreqData : frequency content of the source
% K : number of microphones
% rho : radius of UCA
% L : number of sources to find
% N : length of time domain signal

    estDOA = zeros(L,1); estTOA = zeros(L,1); estBeta = zeros(L,1);
    
    for numSrcs = 1 : L
    
        j0 = numSrcs;
        itErr = 1;
        while ( itErr > 1e-3 )
        
            for j = j0:numSrcs
            
                tmpMicData = micTimeData;
            
                indVec = 1:numSrcs;
                indVec = indVec(indVec ~= j);
            
                for p = indVec
                
                    y_p = genDelayData(srcTimeData,estDOA(p),estTOA(p),K, rho);
                    tmpMicData = tmpMicData - estBeta(p)*y_p;
                
                end
            
                tmpMicFreqData = getFreqMicData(tmpMicData, N, K);
            
                [estDOA(j), estTOA(j), estBeta(j)] = estML_TOA_DOA(tmpMicFreqData,...
                    srcFreqData, K, rho, estDOA(indVec), estTOA(indVec));
%                 [estDOA(j), estTOA(j), estBeta(j)] = estML_TOA_DOA(tmpMicFreqData,...
%                     srcFreqData, K, rho);
                % Bounding beta
%                 estBeta(j) = min(1, estBeta(j));
                estBeta = min(ones(numSrcs,1),updateBeta(micTimeData, srcTimeData, estDOA,...
                    estTOA, N, K, numSrcs, rho));
            end
        
            j0 = 1;
        
            currCost = getCost(micTimeData, srcTimeData,  estDOA,...
                estTOA, estBeta, numSrcs, K, rho);
        
            if (numSrcs == 1)
                prevCost = currCost;
            end
        
%             itErr = abs( currCost - prevCost );
            itErr = 0;
        
            prevCost = currCost;
        
        end
 
    end
end