function [estPos, e, xGrid, yGrid] = extensiveSearch(micTimeData,srcTimeData,...
    micPos,fs, xGrid, yGrid)

    K = size(micPos,1);
    minSeparationOfSources = 1e-3;              % time [s]
    samples2Sum = minSeparationOfSources*fs;    % samples
    
    for ix = 1:length(xGrid)
        for iy = 1:length(yGrid)
            srcTestPos = [xGrid(ix) yGrid(iy) 1.26];
            [~,~,delay] = getMic2SrcParams(micPos',srcTestPos');
            resid = 0;
            
            for k = 1:K
                delaySample = round(delay(k)*fs);
                if mod(delaySample,2) == 0
                    delaySample = delaySample + 1;
                end
                extSignal = [srcTimeData(1:samples2Sum); zeros(delaySample,1)];
                dSignal = nfGenDelayData(extSignal, micPos(k,1:2)', ...
                    srcTestPos(1:2)',fs);
                residTmp = abs(micTimeData(k,1:delaySample+samples2Sum)...
                    - dSignal');
                resid = resid + norm(residTmp);
            end
            
            e(ix,iy) = resid;
        end
    end
    
    [~,maxIndx] = findMaximizer(detrend(e).*detrend(e')');
    estPos = [xGrid(maxIndx(1)) yGrid(maxIndx(2))];

end