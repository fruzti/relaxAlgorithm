% Test of RELAX algorithm estimating L sources, with L known a-priori
clear all, close all

K = 6;
rho = 0.06;
L = 3;

% Creates K virtual signals coming from L random waves impinging in a UCA
% with radius rho.
[micFreqData, srcData, trueDOA, trueTOA, ...
    micTimeData, srcTimeData] = genTstMicData(K, rho,L);

N = length(srcTimeData);

% RELAX : BEGIN

estDOA = zeros(L,1); estTOA = zeros(L,1); estBeta = zeros(L,1);

for numSrcs = 1 : L
    
    j0 = numSrcs;
    
    while ( itErr > 10e-3 )
        
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
                srcData, K, rho);
        end
        
        j0 = 1;
        
        currCost = getCost(micTimeData, estDOA, estTOA, estBeta, numSrcs);
        
        itErr = abs( currCost - prevCost );
        
        prevCost = currCost;
        
    end
 
end