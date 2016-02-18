function [estDOA, estTOA, estBeta, J] = estML_TOA_DOA(micFreqData,srcData,K,rho)
% estDOA : mle of direction-of-arrival
% estTOA : mle of time-of-arrival
% estBeta : mle of amplitude of signal
% micFreqData : fft coefficients of receivers data
% srcData     : fft coefficients of source
% K     : number of microphones
% rho   : radius of UCA

    % number of harmonics
    l = floor(length(srcData)/2);
    
    % uniform grid for cost function
    Theta = deg2rad(0:359); Eta = deg2rad(0:359);
    J = zeros(length(Theta),length(Eta));
    % init maxJ
    maxJ = 0; 
    for i = 1:length(Theta)
        for j = 1:length(Eta)
            J(i,j) = evalCost(micFreqData, Theta(i), srcData, Eta(j), K, l, rho);
            if J(i,j) > maxJ
                maxJ = J(i,j); estDOA = i; estTOA = j;
            end
        end
    end
    % mle of attenuation
    Ps = (2*l + 1) * (srcData'*srcData);
    estBeta = sqrt(maxJ)/(K*Ps);

end