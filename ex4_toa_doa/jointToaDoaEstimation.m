function [estimatedXi, estimatedDoa] = jointToaDoaEstimation(...
        dataMatrix, sourceSignal, rho, nGridPoints, tolerance)
    [nData, nSensors] = size(dataMatrix);
    L = floor(nData/2);
    dftIndices = (-L:L)';
    % compute the symmetric DFT of the source and sensor data
    dftSourceData = dft2SymmetricDft(fft(sourceSignal));
    dftSensorDataMatrix = dft2SymmetricDft(fft(dataMatrix));
    dftCrossCorrelation = (conj(dftSourceData)*ones(1,nSensors)).*...
        dftSensorDataMatrix;
    % set up the DOA and range grid
    nRanges = nGridPoints(1);
    rangeGrid = linspace(0,2*pi,nRanges);
    nDoas = nGridPoints(2);
    doaGrid = linspace(0,2*pi,nDoas);
    % compute the grid cost function
    costFunctionGrid = nan(nRanges,nDoas);
    for iDoa = 1:nDoas
        doaDelayMatrix = exp(-1i*rho*dftIndices*...
            cos(doaGrid(iDoa)-2*pi*(0:nSensors-1)/nSensors));
        modulatedDftCrossCorrelation = ...
            sum(doaDelayMatrix.*dftCrossCorrelation,2);
        costFunctionGrid(:,iDoa) = ...
            abs(myIfft(symmetricDft2Dft(modulatedDftCrossCorrelation, ...
            nData), nRanges)).^2;
    end
    % find the maximum on the grid
    [rangeIndex, doaIndex] = ...
        find(max(max(costFunctionGrid)) == costFunctionGrid);
    coarseRangeEstimate = rangeGrid(rangeIndex(1));
    coarseDoaEstimate = doaGrid(doaIndex(1));
    % refine the DOA estimate
    doaLimits = coarseDoaEstimate+diff(doaGrid(1:2))*[-1;1];
    objectiveFunction = @(doa) -jointToaDoaCostFunction(...
        dftCrossCorrelation, coarseRangeEstimate, doa, L, rho);
        [lowerDoa, upperDoa] = fibonacciSearch(objectiveFunction, ...
            doaLimits(1),doaLimits(2), tolerance(2));
    estimatedDoa = (lowerDoa + upperDoa)/2;
    % refine the range estimate
    rangeLimits = coarseRangeEstimate+diff(rangeGrid(1:2))*[-1;1];
    objectiveFunction = @(xi) -jointToaDoaCostFunction(...
        dftCrossCorrelation, xi, estimatedDoa, L, rho);
        [lowerRange, upperRange] = fibonacciSearch(objectiveFunction, ...
            rangeLimits(1),rangeLimits(2), tolerance(1));
    estimatedXi = (lowerRange + upperRange)/2;
end

