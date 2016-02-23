run experimentalSetup.m

% Simulation Parameters
zPos = 1.26;
srcPos = [ladLoud' zPos];
srcPos(:,1) = srcPos(:,1) + 4;
micPos = [arrayElems; ones(1,size(arrayElems,2))*zPos]';
micPos(:,1) = micPos(:,1) + 4;

arrayPos = mean(micPos(7:12,:));

clearvars -except arrayPos srcPos micPos

N = 81;
srcTimeData = randn(N,1);

fs = 20e3;
K = size(micPos,1);
srcPos2 = [1.5; 0.1; 1.26]';
srcPos3 = [2.6; 0.25; 1.26]';

trueValues = [srcPos(1:2); srcPos2(1:2); srcPos3(1:2)];
%%
% Generation of Signals
[micTimeData, micFreqData, srcFreqData] = nfGenTstData(srcTimeData,...
    micPos, trueValues, fs);
L = size(trueValues,1);

warning off
[estX,estY,estBeta,J] = nfSequentialMLE_TOA_DOA(micTimeData,srcTimeData,...
    srcFreqData,K,micPos,L,N,fs);

xGrid = 0:0.01:4;
yGrid = 0:0.01:0.5;

imagesc(xGrid,yGrid,J{1}'), set(gca,'YDir','normal')
title('Cost Function for Source Location')
xlabel('X-Axis'), ylabel('Y-Axis')
disp('Estimate')
estimates = [estX estY];
[~,or] = sort(estimates(:,1));
disp(estimates(or,:))
disp('True')
[~,or] = sort(trueValues(:,1));
disp(trueValues(or,:))