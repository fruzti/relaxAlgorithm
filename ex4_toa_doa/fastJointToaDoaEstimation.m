function [estimatedXi, estimatedDoa] = fastJointToaDoaEstimation(...
        dataMatrix, sourceSignal, rho, nGridPoints, tolerance)
    [nData, nSensors] = size(dataMatrix);
    L = floor(nData/2);
    dftIndices = (-L:L)';
    % compute the symmetric DFT of the source and sensor data
    dftSourceData = dft2SymmetricDft(fft(sourceSignal));
    dftSensorDataMatrix = dft2SymmetricDft(fft(dataMatrix));
    dftCrossCorrelation = (conj(dftSourceData)*ones(1,nSensors)).*...
        dftSensorDataMatrix;
        
    % we first estimate the DOA
    nDoas = nGridPoints(2);
    doaGrid = linspace(0,2*pi,nDoas);    
    costFunctionGridDoa = nan(nDoas,1);
    for iDoa = 1:nDoas
        doaDelayMatrix = exp(-1i*rho*dftIndices*...
            cos(doaGrid(iDoa)-2*pi*(0:nSensors-1)/nSensors));
        modulatedDftSensorDataMatrix = ...
            sum(doaDelayMatrix.*dftSensorDataMatrix,2);
        costFunctionGridDoa(iDoa) = ...
            real(modulatedDftSensorDataMatrix'*...
            modulatedDftSensorDataMatrix);
    end
    [~, doaIdx] = max(costFunctionGridDoa);
    coarseDoaEstimate = doaGrid(doaIdx(1));
    % refine the doa estimate
    doaLimits = coarseDoaEstimate+diff(doaGrid(1:2))*[-1;1];
    objectiveFunction = @(doa) -doaCostFunction(...
        dftSensorDataMatrix, doa, L, rho);
    [lowerDoa, upperDoa] = fibonacciSearch(objectiveFunction, ...
        doaLimits(1),doaLimits(2), tolerance(2));
    estimatedDoa = (lowerDoa + upperDoa)/2;
    
    % we then find the range
    nRanges = nGridPoints(1);
    rangeGrid = linspace(0,2*pi,nRanges);
    doaDelayMatrix = exp(-1i*rho*dftIndices*...
        cos(estimatedDoa-2*pi*(0:nSensors-1)/nSensors));
    modulatedDftCrossCorrelation = ...
        sum(doaDelayMatrix.*dftCrossCorrelation,2);
    costFunctionGridRange = ...
        abs(myIfft(symmetricDft2Dft(modulatedDftCrossCorrelation, ...
        nData), nRanges)).^2;
    % find the maximum on the grid
    [~,rangeIdx] = max(costFunctionGridRange);
    coarseRangeEstimate = rangeGrid(rangeIdx(1));
    % refine the range estimate
    rangeLimits = coarseRangeEstimate+diff(rangeGrid(1:2))*[-1;1];
    objectiveFunction = @(xi) -jointToaDoaCostFunction(...
        dftCrossCorrelation, xi, estimatedDoa, L, rho);
    [lowerRange, upperRange] = fibonacciSearch(objectiveFunction, ...
        rangeLimits(1),rangeLimits(2), tolerance(1));
    estimatedXi = (lowerRange + upperRange)/2;
end

