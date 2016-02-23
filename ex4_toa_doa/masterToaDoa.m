clear
clc
close all

% add path to fast ML method
addpath ../../fast_f0_estimation/complex/matlab/
% Setup
fastAlgorithmIsIncluded = false;
samplingFreq = 8000;
propagationSpeed = 343;
nData = 100;
nMc = 250;
sensorRadius = 0.1;
rho = 2*pi*sensorRadius*samplingFreq/(nData*propagationSpeed);
nSensors = 3;
gain = 0.9;

% the source signal
sourceSignal = randn(nData,1);

% noise variances
snrList = (-10:1:10)';
nSnrs = length(snrList);
noiseVariances = (sourceSignal'*sourceSignal)/nData*...
    10.^(-snrList/10);

% compute the CRLB for the freqeuncy (in radians/sample)
L = floor(nData/2);
Z = exp(1i*2*pi*(0:nData-1)'*(-L:L)/nData);
alpha = (Z'*sourceSignal)/nData;
weightedAlphas = real(alpha'*(alpha.*((-L:L)'.^2)));
crlbXi = noiseVariances/(nSensors*gain^2*nData*weightedAlphas);
crlbRange = crlbXi*(nData*propagationSpeed/(2*pi*samplingFreq))^2;
crlbDoa = crlbXi*2/rho^2;
crlbDoaDegrees = crlbDoa*(180/pi)^2;

% accuracy of the grid+refinement
refSnr = 15; % dB
accuracyXi = sqrt(10^(-refSnr/10)/(nSensors*gain^2*nData*weightedAlphas));
accuracyDoa = accuracyXi*sqrt(2)/rho;
dbGainList = [10, 3, 1, 0.1, 0.01];
% dbGainList = [0.1];
gainList = 10.^(dbGainList/10);
gridSizeListXi = [sqrt(((gainList-1)./gainList)*...
    (sourceSignal'*sourceSignal/(nData*weightedAlphas))), accuracyXi];
gridSizeListDoa = [2*sqrt(((gainList-1)./gainList)*...
    (sourceSignal'*sourceSignal/(nData*rho^2*weightedAlphas))), accuracyDoa];
nGridSizes = length(gainList)+1;

%% run the Monte Carlo Simulation
estimationErrorXi = nan(nMc, nSnrs, nGridSizes);
estimationErrorDoa = nan(nMc, nSnrs, nGridSizes);
computationTime = nan(nMc, nSnrs, nGridSizes);
if fastAlgorithmIsIncluded
    estimationErrorXiFast = nan(nMc, nSnrs, nGridSizes);
    estimationErrorDoaFast = nan(nMc, nSnrs, nGridSizes);
    computationTimeFast = nan(nMc, nSnrs, nGridSizes);
end
nGridPoints = nan(2, nGridSizes);
for iSnr = 1:nSnrs
    for jMc = 1:nMc
        [snrList(iSnr), jMc]
        expDoa = 2*pi*rand(1); % radians
        expXi = 2*pi*rand(1); % radians
        delays = nData*(expXi-...
            rho*cos(expDoa-2*pi*(0:nSensors-1)/nSensors))/(2*pi);
        % generate sensor data
        cleanDataMatrix = nan(nData,nSensors);
        for iSensor = 1:nSensors
            cleanDataMatrix(:,iSensor) = ...
                gain*delayPeriodicSignal(sourceSignal,delays(iSensor));
        end
        noise = sqrt(noiseVariances(iSnr))*randn(nData,nSensors);
        dataMatrix = cleanDataMatrix+noise;
        for kGridSize = 1:nGridSizes
            nGridPoints(1,kGridSize) = ceil(2*pi/gridSizeListXi(kGridSize));
            nGridPoints(2,kGridSize) = ceil(2*pi/gridSizeListDoa(kGridSize));
            % compute the estimates (standard)
            tic
            [estimatedXi, estimatedDoa] = jointToaDoaEstimation(...
                dataMatrix, sourceSignal, rho, ...
                nGridPoints(:,kGridSize),[accuracyXi, accuracyDoa]);
            computationTime(jMc,iSnr,kGridSize) = toc;
            estimationErrorXi(jMc,iSnr,kGridSize) = ...
                findSmallestAngularError(expXi, estimatedXi);
            estimationErrorDoa(jMc,iSnr,kGridSize) = ...
                findSmallestAngularError(expDoa, estimatedDoa);
            % compute the estimates (fast)
            if fastAlgorithmIsIncluded
                tic
                [estimatedXi, estimatedDoa] = fastJointToaDoaEstimation(...
                    dataMatrix, sourceSignal, rho, ...
                    nGridPoints(:,kGridSize),[accuracyXi, accuracyDoa]);
                computationTimeFast(jMc,iSnr,kGridSize) = toc;
                estimationErrorXiFast(jMc,iSnr,kGridSize) = ...
                    findSmallestAngularError(expXi, estimatedXi);
                estimationErrorDoaFast(jMc,iSnr,kGridSize) = ...
                    findSmallestAngularError(expDoa, estimatedDoa);
            end
        end
    end
end

%% Compute the average error and computation time
rmseXi = nan(nSnrs, nGridSizes);
rmseDoa = nan(nSnrs, nGridSizes);
minComputationTime = nan(nSnrs, nGridSizes);
for kGridSize = 1:nGridSizes
    rmseXi(:,kGridSize) = sqrt(mean(estimationErrorXi(:,:,kGridSize).^2))';
    rmseDoa(:,kGridSize) = sqrt(mean(estimationErrorDoa(:,:,kGridSize).^2))';
    minComputationTime(:,kGridSize) = min(computationTime(:,:,kGridSize))';
end
avgComputationTime = mean(minComputationTime);

if fastAlgorithmIsIncluded
    rmseXiFast = nan(nSnrs, nGridSizes);
    rmseDoaFast = nan(nSnrs, nGridSizes);
    minComputationTimeFast = nan(nSnrs, nGridSizes);
    for kGridSize = 1:nGridSizes
        rmseXiFast(:,kGridSize) = sqrt(mean(estimationErrorXiFast(:,:,kGridSize).^2))';
        rmseDoaFast(:,kGridSize) = sqrt(mean(estimationErrorDoaFast(:,:,kGridSize).^2))';
        minComputationTimeFast(:,kGridSize) = min(computationTimeFast(:,:,kGridSize))';
    end
    avgComputationTimeFast = mean(minComputationTimeFast);
end



%% plot the results

figure(1)
if fastAlgorithmIsIncluded
    semilogy(snrList,[sqrt(crlbXi),rmseXi, rmseXiFast])
    legend(['CRLB';cellstr(num2str(gainList', 'g = %-d'));'No refinement';...
        cellstr(num2str(gainList', 'fast g = %-d'));'Fast no refinement'])
else
    semilogy(snrList,[sqrt(crlbXi),rmseXi])
    legend(['CRLB';cellstr(num2str(gainList', 'g = %-d'));'No refinement'])
end
title('TOA')
xlabel('snr [dB]')
ylabel('RMSE [cycles/sample]')

figure(2)
if fastAlgorithmIsIncluded
    semilogy(snrList,[sqrt(crlbDoa),rmseDoa, rmseDoaFast])
    legend(['CRLB';cellstr(num2str(gainList', 'g = %-d'));'No refinement';...
        cellstr(num2str(gainList', 'fast g = %-d'));'Fast no refinement'])
else
    semilogy(snrList,[sqrt(crlbDoa),rmseDoa])
    legend(['CRLB';cellstr(num2str(gainList', 'g = %-d'));'No refinement'])
end
title('DOA')
xlabel('snr [dB]')
ylabel('RMSE [cycles/sample]')

%% store the data
comment = ['The RMSE vs SNR for various grid sizes. The gains are ', num2str(gainList, 'g = %-d ')];
rmseData = [snrList, sqrt(crlbXi), rmseXi];
labels = 'snr crlb gridA gridB gridC gridD gridE fftOnly';
mtxToTblPgfplots('toa',comment,labels,rmseData)

comment = ['The RMSE vs SNR for various grid sizes. The gains are ', num2str(gainList, 'g = %-d ')];
rmseData = [snrList, sqrt(crlbDoa), rmseDoa];
labels = 'snr crlb gridA gridB gridC gridD gridE fftOnly';
mtxToTblPgfplots('doa',comment,labels,rmseData)