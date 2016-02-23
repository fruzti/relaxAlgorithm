function [estX, estY, estBeta, J] = nfEstML_TOA_DOA(micFreqData,...
    srcFreqData,micPos,fs)
% function [estX, estY, estBeta, J] = nfEstML_TOA_DOA(micFreqData,...
%     srcFreqData,micPos,fs)
% ---------------------------------------------------------------------
% micFreqData : fft coefficients of the microphones signal
% srcFreqData : fft coefficients of the source signal
% micPos : positions of the microphones
% fs : sampling frequency

    % number of samples
    N = length(micFreqData);
    % number of harmonics
    l = floor(N/2);
    % number of microphones
    K = size(micPos,1);
    
%     [xGrid, yGrid] = meshgrid(0:0.01:4, 0:0.01:0.5);
    xGrid = 0:0.01:4;
    yGrid = 0:0.01:0.5;
    
    J = zeros(length(xGrid),length(yGrid));
    maxJ = 0;
    for i = 1:length(xGrid)
        for j = 1:length(yGrid)
            srcPos = [xGrid(i); yGrid(j)];
            J(i,j) = nfEvalCost(micFreqData,srcFreqData, micPos',...
                srcPos,l,N,fs);
            if J(i,j) > maxJ
                maxJ = J(i,j); estX = xGrid(i); estY = yGrid(j);
            end
        end
    end
    
    Ps = (2*l + 1) * (srcFreqData'*srcFreqData);
    estBeta = sqrt(maxJ)/(K*Ps);

end


