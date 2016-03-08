function J = nfEstML(micFreqData,...
    srcFreqData,micPos,fs, xPoints, yPoints)
% function J = nfEstML_TOA_DOA(micFreqData,...
%     srcFreqData,micPos,fs)
% ---------------------------------------------------------------------
% Modifed Version for selected grid (There is no need
% for the estimate in this step
% ---------------------------------------------------------------------
% micFreqData : fft coefficients of the microphones signal
% srcFreqData : fft coefficients of the source signal
% micPos : positions of the microphones
% fs : sampling frequency
% xPoints : x-axis points to evaluate
% yPoints : y-axis points to evaluate

    % number of samples
    N = length(micFreqData);
    % number of harmonics
    l = floor(N/2);
    % number of microphones
    K = size(micPos,1);
    
    % Cost function
    J = zeros(length(xPoints),1);
    
%     maxJ = 0;
    for i = 1:length(xPoints)
        srcPos = [xPoints(i); yPoints(i)];
        J(i) = nfEvalCost(micFreqData,srcFreqData, micPos',...
            srcPos,l,N,fs);
%         if J(i) > maxJ
%             maxJ = J(i); estX = xGrid(i); estY = yGrid(i);
%         end
    end
    
%     myCost = @(x) -nfEvalCost(micFreqData,srcFreqData, micPos',...
%         x,l,N,fs);
%     x0 = [estX; estY];
%     options = optimset('Display', 'off') ;
%     xStar = fminunc(myCost,x0,options);
%     estX = xStar(1); estY = xStar(2);
%     
%     Ps = (2*l + 1) * (srcFreqData'*srcFreqData);
%     estBeta = sqrt(maxJ)/(K*Ps);

end



