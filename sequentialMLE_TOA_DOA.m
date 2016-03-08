function [estDOA, estTOA, estBeta, J] = sequentialMLE_TOA_DOA(micTimeData,...
    srcTimeData, srcFreqData, K, p,L,N)
% function [estDOA, estTOA, estBeta] = sequentialMLE_TOA_DOA(micTimeData,...
%     srcTimeData, srcFreqData, K, p,L,N)
% ------------------------------------------------------------------------
% estDOA : estimated DOAs
% estTOA : estimated TOAs
% estBeta : estimated Betas
% micTimeData : time-domain signal of the receivers
% srcTimeData : time-domain signal from the transmitter
% srcFreqData : frequency content of the source
% K : number of microphones
% p : frequency-radius equivalence
% L : number of sources to find
% N : length of time domain signal
% J : cost function per iteration

    estDOA = zeros(L,1); estTOA = zeros(L,1); estBeta = zeros(L,1);
    options = optimset('Display', 'off') ;
    J = cell(L,1);

    for numSrcs = 1 : L
    
        j0 = numSrcs;
        itErr = 1;
        
        while ( itErr > 1e-3 )
        
            for j = j0:numSrcs
            
                tmpMicData = micTimeData;
            
                indVec = 1:numSrcs;
                indVec = indVec(indVec ~= j);
                
                for indx = indVec
                
                    y_p = genDelayData(srcTimeData,estDOA(indx),estTOA(indx),K, p);
                    tmpMicData = tmpMicData - estBeta(indx)*y_p;
                
                end
            
                tmpMicFreqData = getFreqMicData(tmpMicData, N, K);
            
                [estDOA(j), estTOA(j), ~, J{j}] = estML_TOA_DOA(tmpMicFreqData,...
                    srcFreqData, K, p, estDOA(indVec), estTOA(indVec));
%                 [estDOA(j), estTOA(j), ~, J{j}] = estML_TOA_DOA(tmpMicFreqData,...
%                     srcFreqData, K, p);
                
                % Refining the grid [Line Search]
                myCost = @(x) -evalCost(tmpMicFreqData,x(1),...
                    srcFreqData,x(2),K,floor(size(tmpMicFreqData,1)/2),p);
                x0 = [estDOA(j); estTOA(j)];
                xStar = fminunc(myCost,x0,options);
                
                % Updating with finer grid
                estDOA(j) = xStar(1); estTOA(j) = xStar(2);
                
                % Bounding beta to unit
                estBeta = min(ones(numSrcs,1),updateBeta(micTimeData, srcTimeData, estDOA,...
                    estTOA, N, K, numSrcs, p));
            end
        
            j0 = 1;
        
            currCost = getCost(micTimeData, srcTimeData,  estDOA,...
                estTOA, estBeta, numSrcs, K, p);
        
            if (numSrcs == 1)
                prevCost = currCost;
            end
        
            % RELAX Mode
%             itErr = abs( currCost - prevCost );

            % Sequential Mode
            itErr = 0;
        
            prevCost = currCost;
        
        end
 
    end
end