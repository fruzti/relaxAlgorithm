function value = jointToaDoaCostFunction(dftCrossCorrelation, xi, doa, ...
        L, rho)
    nSensors = size(dftCrossCorrelation,2);
    angularDelays = xi-rho*cos(doa-2*pi*(0:nSensors-1)/nSensors);
    delayVectors = exp(1i*(-L:L)'*angularDelays);
    value = abs(sum(sum(delayVectors.*dftCrossCorrelation,1),2))^2;
end

