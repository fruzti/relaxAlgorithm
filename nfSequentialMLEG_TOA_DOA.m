function [estX, estY, estBeta, J] = nfSequentialMLEG_TOA_DOA(micTimeData,...
    srcTimeData, srcFreqData, K, micPos,L,N,fs,mask)
% function [estX, estY, estBeta, J] = nfSequentialMLE_TOA_DOA(micTimeData,...
%     srcTimeData, srcFreqData, K, micPos,L,N,fs)
% ------------------------------------------------------------------------
% estX : estimated x-axis positions
% estY : estimated y-axis positions
% estBeta : estimated Betas
% J : cost function per iteration
% micTimeData : time-domain signal of the receivers
% srcTimeData : time-domain signal from the transmitter
% srcFreqData : frequency content of the source
% K : number of microphones
% micPos : microphone positions
% L : number of sources to find
% N : length of time domain signal
% fs : sampling frequency

    estX = zeros(L,1); estY = zeros(L,1); estBeta = zeros(L,1);
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
                    tmpSrcPos = [estX(indx); estY(indx)];
                    y_p = nfGenDelayData(srcTimeData,micPos(:,1:2)',tmpSrcPos,fs);
                    tmpMicData = tmpMicData - estBeta(indx)*y_p;
                
                end
            
                tmpMicFreqData = getFreqMicData(tmpMicData, N, K);
            
                
                [estX(j), estY(j), ~, J{j}] = nfEstMLG_TOA_DOA(tmpMicFreqData,...
                    srcFreqData, micPos(:,1:2), fs,mask);

                % Bounding beta to unit
                estBeta = min(ones(numSrcs,1),...
                    nfUpdateBeta(micTimeData, srcTimeData, micPos(:,1:2)',...
                    estX, estY, N, K, numSrcs,fs));
            end
        
            j0 = 1;
        
            currCost = nfGetCost(micTimeData,srcTimeData,micPos(:,1:2)',estX,...
                estY,estBeta,numSrcs,fs);
        
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