function value = doaCostFunction(dftSensorDataMatrix, doa, ...
        L, rho)
    nSensors = size(dftSensorDataMatrix,2);
    angularDelays = -rho*cos(doa-2*pi*(0:nSensors-1)/nSensors);
    delayVectors = exp(1i*(-L:L)'*angularDelays);
    value = sum(abs(sum(delayVectors.*dftSensorDataMatrix,2)).^2);
end

