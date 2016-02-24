run experimentalSetup.m

% Simulation Parameters
zPos = 1.26;
srcPos = [ladLoud' zPos];
srcPos(:,1) = srcPos(:,1) + 4;
micPos = [arrayElems; ones(1,size(arrayElems,2))*zPos]';
micPos(:,1) = micPos(:,1) + 4;

arrayPos = mean(micPos(7:12,:));

clearvars -except arrayPos srcPos micPos
%%
% N = 1921;
N = 101;
% N = 501;
% srcTimeData = randn(N,1);
srcTimeData = mls(7,1);
% srcTimeData = repmat([srcTimeData; zeros(30,1)],3,1);
N = length(srcTimeData);


fs = 10e3;
K = size(micPos,1);
% srcPos2 = [1.5; 0.1; 1.26]';
% srcPos3 = [2.6; 0.25; 1.26]';

vrtSrcPos = genSrcsFromWalls(srcPos(1:2), [4 6]);
vrtSrcPos = vrtSrcPos(:,[1 3:4]);
trueValues = [srcPos(1:2); vrtSrcPos'];
% trueValues = [srcPos(1:2); srcPos2(1:2); srcPos3(1:2)];
%%
srcTimeData = mls(7,4); N1 = 128;
% srcTimeData = randn(201,1);
fileNameRIRs = 'ISM_RIRs.mat';
AuData = ISM_AudioData(fileNameRIRs, srcTimeData);
srcTimeData = srcTimeData(N1:end);
micTimeData = AuData(N1:end,1:18); N = size(micTimeData,1); fs = 10e3;
K = size(micPos,1);
%%
fs = 5e3;
micTimeData = resample(micTimeData,fs,10e3);
srcTimeData = resample(srcTimeData,fs, 10e3); N = size(micTimeData,1);
%%
micFreqData = getFreqMicData(micTimeData, N, K);
srcFreqData = applyFFT(srcTimeData, N);
L = 4;
%%
% Generation of Signals
% [micTimeData, micFreqData, srcFreqData] = nfGenTstData(srcTimeData,...
%     micPos, trueValues, fs);
% L = size(trueValues,1);

warning off
[estX,estY,estBeta,J] = nfSequentialMLE_TOA_DOA(micTimeData,srcTimeData,...
    srcFreqData,K,micPos,L,N,fs);

% xGrid = 0:0.01:4;
% yGrid = 0:0.01:0.5;
% xGrid = -1:0.1:8;
% yGrid = 0:0.1:12;
xGrid = -3:0.1:8;
yGrid = -3:0.1:12;

figure, imagesc(xGrid,yGrid,J{1}'), set(gca,'YDir','normal')
title('Cost Function for Source Location')
xlabel('X-Axis'), ylabel('Y-Axis')
disp('Estimate')
estimates = [estX estY];
[~,or] = sort(estimates(:,1));
disp(estimates(or,:))
disp('True')
[~,or] = sort(trueValues(:,1));
disp(trueValues(or,:))