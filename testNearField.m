run experimentalSetup.m

% Simulation Parameters
zPos = 1.26;
srcPos = [ladLoud' zPos];
srcPos(:,1) = srcPos(:,1) + 4;
micPos = [arrayElems; ones(1,size(arrayElems,2))*zPos]';
micPos(:,1) = micPos(:,1) + 4;

arrayPos = mean(micPos(7:12,:));

clearvars -except arrayPos srcPos micPos

N = 201;
srcTimeData = randn(N,1);

fs = 10e3;
K = size(micPos,1);

micTimeData = nfGenDelayData(srcTimeData, micPos', srcPos', fs);

micFreqData = getFreqMicData(micTimeData, N, K);
srcFreqData = applyFFT(srcTimeData, N);

[estX, estY, estBeta, J] = nfEstML_TOA_DOA(micFreqData, srcFreqData,...
    micPos(:,1:2), fs);
xGrid = 0:0.01:4;
yGrid = 0:0.01:0.5;
imagesc(xGrid,yGrid,J'), set(gca,'YDir','normal')
title('Cost Function for Source Location')
xlabel('X-Axis'), ylabel('Y-Axis')
disp('Estimate')
disp([estX estY])
disp('True')
disp(srcPos(1:2))