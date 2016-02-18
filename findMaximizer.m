function [maxValue, maxIndx] = findMaximizer(J)
% maxValue : maximum of cost function
% maxIndx : index of the maximum of cost function
% J : cost function

    dims = 2;   % hardcoded dimensionality (two parameters assumed)
    maxIndx = zeros(dims,1);
    
    [maxTmp, maxRows] = max(J);
    [maxValue, maxIndx(2)] = max(maxTmp);
    maxIndx(1) = maxRows(maxIndx(2));

end